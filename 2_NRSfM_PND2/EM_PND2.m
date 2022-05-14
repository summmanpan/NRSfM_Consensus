function [X, s, R, M, C, var_D] = EM_PND2(s, R, X, M, D, W, Gth, Mth, iter_max)
% function [X, s, R, M, C, var_D] = EM_PND2(s, R, X, M, D, W, Gth, Mth, iter_max)
%
% Fit PND to a given data using EM algorithm
%
% Common:
%     s: Scales       (1 x n)
%     R: Rotations    (3 x 3 x n)
%     M: Mean shape   (3 x p)
%     X: Reconstructed shapes (3 x p x n)
%     
% Inputs:
%     D: Observations                                           (3 x p x n)
%     W: Weights - observed (true) or missing (false) (3 x p x n)
%     Gth: Threshold for iteration (cost)
%     Mth: Threshold for iteration (change of mean shape)
%     iter_max: Maximum number of iteration
%
% Outputs:
%     C: Shape covariance     (3p x 3p)
%     var_D: Noise variance
%
% Ref: Minsik Lee, Jungchan Cho, and Songhwai Oh,
% "Procrustean Normal Distribution for Non-Rigid Structure from Motion,"
% IEEE Trans. Pattern Analysis and Machine Intelligence, to appear.
%
% Authors: Minsik Lee (mlee.paper@gmail.com)
% Last update: 2016-09-08
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


% initial values and parameters
C_init = 1e-4;
var_D = 1e-4;
alpha = 1e-4;
v_comp = 2;

if nargin < 7 || isempty(th)
    Gth = 1e-2;
end
if nargin < 8 || isempty(th)
    Mth = eps;
end
if nargin < 9 || isempty(iter_max)
    iter_max = 1e3;
end

[k, p, nSample] = size(D);

pq = p-1;
tdim = k*p;
pdim = k*pq;
ndim = k+1;
cdim = pdim-ndim;
[Qt, ~] = qr(ones(p, 1));
Qt = Qt(:, 2:end);

WW = zeros(pq, pq, k, k, nSample);
for i=1:nSample
    for j=1:k
        tW = W(:, :, i);
        sW = sum(tW(j, :));
        if sW == p
            WW(:, :, j, j, i) = eye(pq);
        elseif sW > 0
            tQ = Qt(tW(j, :), :);
            t2 = sum(tQ)/sqrt(sW);
            WW(:, :, j, j, i) = tQ'*tQ - t2'*t2;
        end
    end
end
WW = reshape(permute(WW, [3 1 4 2 5]), pdim, pdim, nSample);
dW = sum(W(:))-sum(sum(any(W, 2)));

% Exclude translation bases
M = M*Qt;
X = ipermute(reshape(reshape(permute(X, [3 1 2]), [], p)*Qt, [], k, pq), [3 1 2]);
wX = zeros(k, pq, nSample);
qD = reshape(ipermute(reshape(reshape(permute(D, [1 3 2]), k*nSample, p)*Qt, k, nSample, pq), [1 3 2]), pdim, nSample);

sR = bsxfun(@times, R, reshape(s, 1, 1, []));
sRX = squeeze(sum(bsxfun(@times, reshape(sR, k, k, 1, []), reshape(X, 1, k, pq, [])), 2));

% Caculate initial Q
temp = zeros(k, pq, k);
temp(1, :, 2) = -M(3, :);
temp(1, :, 3) = M(2, :);
temp(2, :, 1) = M(3, :);
temp(2, :, 3) = -M(1, :);
temp(3, :, 1) = -M(2, :);
temp(3, :, 2) = M(1, :);
[Q, ~] = qr([reshape(temp, pdim, k) M(:)]);
Q = Q(:, ndim+1:end);

% Initial shape covariance
vX = Q'*reshape(sRX, [], nSample);
Cp = eye(cdim)*C_init + vX*vX'/nSample;
tC = Q/(Cp/var_D)*Q';
C_X = zeros(pdim, pdim, nSample);


count = 0;
cost = inf;
eG = inf;
eM = inf;
while eG > Gth && eM > Mth && count < iter_max
    pcost = cost;
    pM = M;

    %%% E-Step %%%
    for i=1:nSample
        C_X(:, :, i) = inv(reshape(sR(:, :, i)'*reshape(reshape(sR(:, :, i)'*reshape(tC, k, []), pdim, pdim)', k, []), pdim, pdim) + WW(:, :, i));
        tmp = C_X(:, :, i)*qD(:, i);
        X(:, :, i) = reshape(tmp, k, pq);
        wX(:, :, i) = reshape(WW(:, :, i)*tmp, k, pq);
    end
    C_X = C_X*var_D;
    
    %%% M-Step %%%

    % var_D
    pv = var_D;
    var_D = (norm(qD(:)-wX(:))^2 + WW(:)'*C_X(:))/dW;
    var_D = min(var_D*v_comp, max(pv, var_D));
    var_D = max(var_D, 1e-10);

    % M, s, R, Q, Cp
    C_X = C_X + bsxfun(@times, reshape(X, pdim, 1, []), reshape(X, 1, pdim, []));
    [M, s, R, sR, Q, Cp, G, talpha, tcost] = update_parameters(X, C_X, M, alpha);
    alpha = max(talpha, alpha*(~talpha));
    tC = G*G'*var_D;
    
    %%% Covergence check / calculating the cost value %%%
    cost = (tcost*2 + log(var_D)*dW/nSample)/cdim;
    eG = pcost - cost;
    eM = norm(M(:) - pM(:))^2;
    
    count = count + 1;
    disp(['EM-PND2 ' num2str(count) ' : eM = ' num2str(eM) ' / var_D = ' num2str(var_D) ' / e_cost = ' num2str(eG) ' / cost = ' num2str(cost)]);
end

% Restore translation bases
M = M*Qt';
X = ipermute(reshape(reshape(permute(X, [1 3 2]), [], pq)*Qt', k, nSample, p), [1 3 2]);
C = reshape(permute(reshape(Qt*reshape(permute(reshape(Qt*reshape(permute(reshape(Q*Cp*Q', k, pq, k, pq), [2 3 4 1]), pq, []), p, k, pq, k), [3 4 1 2]), pq, []), p, k, p, k), [2 3 4 1]), tdim, tdim);
C = (C+C')/2;

end

