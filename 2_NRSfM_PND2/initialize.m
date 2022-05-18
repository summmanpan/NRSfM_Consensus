function [D, nD, s, R, X, M] = initialize(D, W, rrth)
% function [D, nD, s, R, X, M] = initialize(D, W, rrth)
%
% Initialize rotations based on orthonormality constraint
%
% Common:
%     D: Observations                                 (3 x p x n)
% 
% Inputs:
%     W: Weights - observed (true) or missing (false) (3 x p x n)
%     rrth: Threshold for bases
%     K_max: Maximum number of K
%
% Outputs:
%     nD: Normalization constant
%     s: Scales                                       (1 x n)
%     R: Rotations                                    (3 x 3 x n)
%     X: Initial shapes                               (3 x p x n)
%     M: Mean shape                                   (3 x p)
%
% Ref: Minsik Lee, Jungchan Cho, and Songhwai Oh,
% "Procrustean Normal Distribution for Non-Rigid Structure from Motion,"
% IEEE Trans. Pattern Analysis and Machine Intelligence, to appear.
%
% Authors: Minsik Lee (mlee.paper@gmail.com)
% Last update: 2016-09-07
% License: GPLv3

%
% Copyright (C) 2016 Minsik Lee
% This file is part of NRSfM_PND2.
%
% NRSfM_PND2 is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% NRSfM_PND2 is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with NRSfM_PND2.  If not, see <http://www.gnu.org/licenses/>.


if nargin < 3 || isempty(rrth)
    rrth = 1e-5;
end

[k, p, nSample] = size(D);
tdim = k*p;

% Remove translation components / normalize
D = pout_trans(D, W); % AQUI USAN W 
nD = sqrt(mse(D(W))*tdim);
D = D/nD;

% Initial scale
s = 1./sqrt(reshape(sum(sum(D.^2, 1), 2), 1, nSample));
fD = bsxfun(@times, permute(D(1:2, :, :), [1 3 2]), s);

% Fill missing entries if needed
if ~all(all(W(1:2, :)))
    T = reshape(permute(~W(1:2, :, :), [1 3 2]), [], p);
    fD = reshape(fill_missing(reshape(fD, [], p), T), 2, nSample, p);
    ts = 1./sqrt(sum(sum(fD.^2), 3));
    s = s.*ts;
    fD = bsxfun(@times, fD, ts);
end

% Calculate rotations
[R, ts] = cal_rotation(reshape(fD, [], p), 1-rrth);
s = s.*ts;
fD = bsxfun(@times, fD, ts);

% Reflection correction
Y = squeeze(sum(bsxfun(@times, reshape(R(:, 1:2, :), 3, 2, 1, nSample), reshape(permute(fD, [1 3 2]), 1, 2, p, nSample)), 2));
[Y, R] = correct_reflection(Y, R);

% Calculate initial shapes
Y = init_shape(Y, R(:, 3, :));

% Transform back to camera coordinate system
X = squeeze(sum(bsxfun(@times, reshape(bsxfun(@rdivide, R, reshape(s, 1, 1, nSample)), k, k, 1, nSample), reshape(Y, k, 1, p, nSample))));

% (Modified) GPA alignment / mean calculation
[~, ts, tR, M] = GPTA(Y);
s = s.*ts;
R = squeeze(sum(bsxfun(@times, reshape(tR, k, k, 1, nSample), reshape(R, 1, k, k, nSample)), 2));

end


function Y = init_shape(X, R)
% Calculate initial shapes based on minimizing the trace of sample covariance
%
% X: Input shapes (R*fD)
% R: The third axes of rotations
% Y: Intial shapes

nSample = size(X, 3);

tR = reshape(R, [], nSample);
tmp = reshape(pout_trans(inner(R, pout_sample(X))), [], nSample);
Z = reshape(tmp*tR'/(tR*tR' - nSample*eye(size(tR, 1)))*tR - tmp, 1, [], nSample);

Y = X + outer(R, Z);

end


function R = outer(X, Y)
% Outer product

    R = bsxfun(@times, X, Y);
end

function R = inner(X, Y)
% Inner product

    R = sum(bsxfun(@times, X, Y));
end


function X = pout_sample(X)
% Eliminate mean component

    X = bsxfun(@minus, X, mean(X, 3));

end

function X = pout_trans(X, W)
% Eliminate translation component

if nargin < 2
    X = bsxfun(@minus, X, mean(X, 2));
else
    X = bsxfun(@minus, X, sum(X.*W, 2)./sum(W, 2)).*W;
    X(isnan(X)) = 0;
end

end


function [X, R] = correct_reflection(X, R)
% Correct reflections signs between different frames
%
% X: 2D shapes (aligned in 3D)
% R: Rotation matrices

nSample = size(X, 3);

pr = true(1, nSample);
r = false(1, nSample);
while any(xor(r, pr))
    pr = r;
    mX = mean(X, 3);
    ind = sum(sum(bsxfun(@times, X, mX))) < 0;
    X(:, :, ind) = -X(:, :, ind);
    r(ind) = ~r(ind);
end
R(:, 1:2, r) = -R(:, 1:2, r);

for i=1:nSample
    if det(R(:, :, i)) < 0
        r(i) = true;
    else
        r(i) = false;
    end
end
R(:, 3, r) = -R(:, 3, r);

end


function [R, s] = cal_rotation(X, th)
% Non-rigid factorization (rotation calculation)
%
% X: (kf x p) data matrix
% th: Threshold for maximum K
% R: Rotations
% s: Scales
%
% This is a modified version of the rotation calculation scheme in
% P. F. U. Gotardo and A. M. Martinez,
% "Non-Rigid Structure from Motion with Complementary Rank-3 Spaces,"
% CVPR 2011, Colorado Springs, Colorado, June 20-25, 2011.

[U, S, ~] = svd(X, 'econ');
ss = diag(S);
ss = cumsum(ss(1:end-1).^2);
maxK = sum(ss/ss(end) < th)+1;
maxK = min(ceil(maxK/3)*3, numel(ss));
F = U(:, 1:maxK);

lcost = inf;
for K = [3:3:size(F, 2)-1 size(F, 2)]
    plcost = lcost;

    [AA1, AA2, off] = cal_A_parts(F(:, 1:K));
    [A, LA] = pre([AA1-AA2; off]);
    if K == 3
        [GG, G] = cal_GG(A, eye(3)/3, 1e-12);
    else
        [GG, G] = cal_GG_m(A, [G zeros(size(G, 1), 3-size(G, 2)); zeros(K-size(G, 1), 3)], 3);
    end

    lcost = norm(A*GG(:))^2*LA;
    
%     disp([num2str(K) ' : ' num2str(plcost-lcost) ' / ' num2str(lcost)]);
    
    if lcost < eps
        break;
    end
end

[R, s] = F2R(F(:, 1:K)*[G zeros(size(G, 1), 3-size(G, 2))]);

end


function [R, s] = F2R(F)
% Project to rotation matrices

    s = zeros(1, size(F, 1)/2);
    R = zeros(3, 3, numel(s));
    R(:, 1:2, :) = reshape(F', 3, 2, []);
    for i=1:numel(s)
        [U, S, V] = svd(R(:, :, i));
        R(:, :, i) = U*V';
        s(i) = 2/trace(S);
    end
    s = s/min(s);

end


function [A, LA] = pre(A)
% Precondition A matrix

    if diff(size(A)) < 0
        [~, A] = qr(A, 0);
    end
    LA = max(eig(A*A'));
    A = A/sqrt(LA);
end

function [GG, G, r] = cal_GG(A, GG, th)
% APG method for finding "rectification" matrix
%
% A: Basis matrix
% GG: outer product of G
% th: Threshold
% G: Rectification matrix
% r: Final rank of G

    if nargin < 3 || isempty(th)
        th = 1e-10;
    end
    cnt = 0;
    vcost = inf;
    AGG = A*GG(:);
    cost = norm(AGG)^2;
    tGG = GG;
    tAGG = AGG;
    t = 1;
    vGG = zeros(size(GG));
    while vcost/cost - 1 > th && vcost > eps
        pGG = GG;
        pAGG = AGG;
        pt = t;

        vGG(:) = tGG(:) - A'*tAGG;
        [GG, G, r] = prj_GG_u(vGG);

        AGG = A*GG(:);
        cost = norm(AGG)^2;
        
        vcost = norm(GG(:)-tGG(:))^2 + cost - norm(tAGG-AGG)^2;

        t = (1+sqrt(1+4*pt^2))/2;
        tGG = GG + (pt-1)/t*(GG-pGG);
        tAGG = AGG + (pt-1)/t*(AGG-pAGG);

        cnt = cnt+1;
%         if mod(cnt, 1e4) == 0
%             disp(['cal_GG ' num2str(cnt) ' : ' num2str(vcost - cost) ' / ' num2str(vcost) ' / ' num2str(cost) ' / ' num2str(norm(pGG(:)-GG(:))^2) ' / ' num2str(r) ' / ' num2str(length(GG))]);
%         end
    end
%     disp(['cal_GG ' num2str(cnt) ' : ' num2str(vcost - cost) ' / ' num2str(vcost) ' / ' num2str(cost) ' / ' num2str(norm(pGG(:)-GG(:))^2) ' / ' num2str(r) ' / ' num2str(length(GG))]);
end

function [GG, G, r] = cal_GG_m(A, G, m, th)
% PG method for finding "rectification" matrix (when K > 3)
%
% A: Basis matrix
% G: Rectification matrix
% GG: outer product of G
% m: Rank constraint
% th: Threshold
% r: Final rank of G

    if nargin < 4 || isempty(th)
        th = 1e-6;
    end
    G = fro_normalize(G(:, 1:m));
    GG = G*G';
    cnt = 0;
    pcost = inf;
    AGG = A*GG(:);
    cost = norm(AGG)^2;
    tID = tic;
    while pcost/cost - 1 > th
        pcost = cost;
        pGG = GG;

        GG(:) = GG(:) - A'*AGG;
        [GG, G, r] = prj_GG_u(GG, m);

        AGG = A*GG(:);
        cost = norm(AGG)^2;

        cnt = cnt+1;
%         if toc(tID) > 1
%             disp(['cal_GG_m ' num2str(cnt) ' : ' num2str(pcost/cost - 1) ' / ' num2str(cost) ' / ' num2str(norm(pGG(:)-GG(:))^2) ' / ' num2str(r) ' / ' num2str(m)]);
%             tID = tic;
%         end
    end
%     disp(['cal_GG_m ' num2str(cnt) ' : ' num2str(pcost/cost - 1) ' / ' num2str(cost) ' / ' num2str(norm(pGG(:)-GG(:))^2) ' / ' num2str(r) ' / ' num2str(m)]);
end


function G = fro_normalize(G)
% Normalize based on Frobenius norm

    G = G/norm(G(:));
end


function [AA1, AA2, off] = cal_A_parts(U)
% Calculate basis matrix for rotation calculation

K = size(U, 2);

A1 = U(1:2:end, :);
A2 = U(2:2:end, :);

AA1 = kron(A1, ones(1, K)).*repmat(A1, 1, K);
AA2 = kron(A2, ones(1, K)).*repmat(A2, 1, K);
off = repmat(reshape(A1, [], 1, K), 1, K).*repmat(A2, [1 1 K]);
off = reshape(off + permute(off, [1 3 2]), size(off, 1), []);

end


function [GG, G, r] = prj_GG_u(GG, m)
% Project to rank-m PSD matrix with unit trace
%
% GG: Input / output matrix
% m: Rank constraint
% G: Decomposed matrix
% r: Final rank of G

    [V, S] = eig(GG);
    if nargin < 2
        [s, ind] = prj_s_u(diag(S));
    else
        [s, ind] = prj_s_u(diag(S), m);
    end
    G = bsxfun(@times, V(:, ind), sqrt(s)');
    GG = G*G';
    r = sum(s > 0);
end


function [s, ind] = prj_s_u(s, m)
% Project to nonnegative vector with m nonzero elements and unit sum
%
% s: Input / output vector (output is sorted)
% m: Maximum number of nonzero elements
% ind: order of elements by sort

    [s, ind] = sort(s, 'descend');
    if nargin == 2
        s = s(1:m);
    end
    dcs = cumsum(s)-s.*(1:numel(s))';
    k = sum(dcs < 1);
    s = s(1:k) - s(k) + (1-dcs(k))/k;
    ind = ind(1:k);
end


function X = fill_missing(Y, T)
% Fill missing entries based on matrix completion
%
% Y: Input
% T: Indicator matrix
% X: Output

isT = 1./(size(T, 2) - sum(T, 2));
isT(isinf(isT)) = 0;

mu = 1e-0;
rho = 1.05;

count = 0;
cost = inf;
tX = Y;
X = zeros(size(Y));
L = zeros(size(Y));
while cost > 1e-10

    [U, S, V] = svd(tX - L, 'econ');
    s = max(diag(S) - 1/mu, 0);
    X = U*diag(s)*V';
    
    temp = (pout_trans(X+L)-Y).*T;
    temp = bsxfun(@plus, temp, sum(temp, 2).*isT);
    M = temp(T);
    
    tX = Y;
    tX(T) = M;
    tX = pout_trans(tX);
    
    Q = X - tX;
    L = L + Q;
    cost = norm(Q(:))^2;

    mu = mu*rho;
    L = L/rho;

    count = count + 1;
%     if mod(count, 1e2) == 0
%         disp(['fill ' num2str(count) ' : ' num2str(sum(s)) ' / ' num2str(cost) ' / ' num2str(mu)]);
%     end
end
% disp(['fill ' num2str(count) ' : ' num2str(sum(s)) ' / ' num2str(cost) ' / ' num2str(mu)]);

end

