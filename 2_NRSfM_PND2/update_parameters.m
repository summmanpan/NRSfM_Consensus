function [m, s, R, sR, Q, Hp, B, alpha, cost] = update_parameters(X, H, m, alpha)
% function [m, s, R, sR, Q, Hp, B, alpha, cost] = update_parameters(X, H, m, alpha)
%
% Update parameters based on gradient descent + line search
%
% Common:
%     m: Mean shape   (3 x p)
%     alpha: Step size
%     
% Inputs:
%     X: Reconstructed shapes           (3 x p x n)
%     H: Convariances of X's (E-step)   (3p x 3p x n)
%
% Outputs:
%     s: Scales                                             (1 x n)
%     R: Rotations                                          (3 x 3 x n)
%     sR: Scaled rotations                                  (3 x 3 x n)
%     Q: Deformation subspace                               (3p x 3p)
%     Hp: Reduced covariance                                ((3p-7) x (3p-7))
%     B: Cholesky decomposition of inverted convariance     (3p x 3p)
%     cost: Cost value
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

[k, p, n] = size(X);
pdim = k*p;
cdim = pdim - k - 1;

if nargin < 4
    alpha = 1e-4;
end
eta = 0.5;
c1 = 1e-4;
c2 = 0.9;

% Calculate gradient
[G, s, R, sR, Q, Hp, B, cost] = find_G(m);
nG = norm(G(:))^2;
c1G = c1*nG;
c2G = c2*nG;

if nG < eps
    return;
end

% Simple line search based on strong Wolfe conditions
maxalpha = alpha;
minalpha = 0;
while true
    tm = m - alpha*G;
    tm = tm/norm(tm(:));
    [tG, ts, tR, tsR, tQ, tHp, tB, tcost] = find_G(tm);
    if tcost > cost-alpha*c1G
        maxalpha = alpha;
        alpha = minalpha + eta*(maxalpha-minalpha);
    else
        pp = G(:)'*tG(:);
        if pp > c2G
            minalpha = alpha;
            if maxalpha == alpha
                maxalpha = alpha/eta;
                alpha = maxalpha;
            else
                alpha = minalpha + eta*(maxalpha-minalpha);
            end
        elseif abs(pp) > c2G
            maxalpha = alpha;
            alpha = minalpha + eta*(maxalpha-minalpha);
        else
            break;
        end
    end
    if maxalpha-minalpha < eps
        alpha = 0;
        break;
    end
end
m = tm;
s = ts;
R = tR;
sR = tsR;
Q = tQ;
Hp = tHp;
B = tB;
cost = tcost;

    function [G, s, R, sR, Q, Hp, B, cost] = find_G(m)
        s = zeros(1, n);
        R = zeros(k, k, n);
        XM = squeeze(sum(bsxfun(@times, reshape(X, k, p, 1, n), reshape(m', 1, p, k)), 2));
        for i=1:n
            [U, S, V] = svd(XM(:, :, i));
            R(:, :, i) = V*U';
            s(i) = 1/trace(S);
        end

        sR = bsxfun(@times, R, reshape(s, 1, 1, []));
        sRX = squeeze(sum(bsxfun(@times, reshape(sR, k, k, 1, []), reshape(X, 1, k, p, [])), 2));
        
        M = zeros(k, p, k);
        M(1, :, 2) = -m(3, :);
        M(1, :, 3) = m(2, :);
        M(2, :, 1) = m(3, :);
        M(2, :, 3) = -m(1, :);
        M(3, :, 1) = -m(2, :);
        M(3, :, 2) = m(1, :);
        M = [m(:) reshape(M, [], k)];
        
        [Q, ~] = qr(M);
        Q = Q(:, size(M, 2)+1:end);
        
        sHc = reshape(sum(bsxfun(@times, sum(bsxfun(@times, reshape(sR, k, k, 1, 1, 1, 1, n), reshape(H, 1, k, p, k, 1, p, n)), 2), permute(sR, [4 5 6 2 1 7 3])), 4), pdim, pdim, n);
        sH = sum(sHc, 3);
        
        Hp = Q'*sH*Q;
        Hp = (Hp+Hp')/(2*n);
        cHp = chol(Hp);
        B = Q/cHp;
        
        G = zeros(size(m));
        BB = B*B';
        
        tmp = reshape(pinv(M)*sH*BB, k+1, k, p);

        G(1, :) = G(1, :) + reshape(- tmp(1, 1, :) - tmp(3, 3, :) + tmp(4, 2, :), 1, []);
        G(2, :) = G(2, :) + reshape(- tmp(1, 2, :) + tmp(2, 3, :) - tmp(4, 1, :), 1, []);
        G(3, :) = G(3, :) + reshape(- tmp(1, 3, :) - tmp(2, 2, :) + tmp(3, 1, :), 1, []);
        
        MX = permute(sum(bsxfun(@times, reshape(sR, k, k, 1, n), reshape(XM, 1, k, k, n)), 2), [3 1 4 2]);
        trMX = MX(1, 1, :) + MX(2, 2, :) + MX(3, 3, :);
        psit = zeros(4, 4, n);
        psit(1, 1, :) = -trMX;
        psit(2, 1, :) = MX(3, 2, :) - MX(2, 3, :);
        psit(3, 1, :) = MX(1, 3, :) - MX(3, 1, :);
        psit(4, 1, :) = MX(2, 1, :) - MX(1, 2, :);   
        psit(1, 2:4, :) = psit(2:4, 1, :);
        psit(2:4, 2:4, :) = -MX;
        psit(2, 2, :) = psit(2, 2, :) + trMX;
        psit(3, 3, :) = psit(3, 3, :) + trMX;
        psit(4, 4, :) = psit(4, 4, :) + trMX;
        
        E = squeeze(sum(bsxfun(@times, reshape(sHc, k, [], 1, n), reshape(reshape(BB, [], p)', 1, [], k)), 2));
        EL = zeros(4, n);
        EL(1, :) = E(1, 1, :) + E(2, 2, :) + E(3, 3, :) - cdim;
        EL(2, :) = E(3, 2, :) - E(2, 3, :);
        EL(3, :) = E(1, 3, :) - E(3, 1, :);
        EL(4, :) = E(2, 1, :) - E(1, 2, :);
        
        O = zeros(k+1, n);
        for i=1:n
            O(:, i) = psit(:, :, i)\EL(:, i);
        end
        
        G = G + reshape(reshape(sRX, [], n)*O(1, :)', k, p);
        G(1, :) = G(1, :) + (squeeze(sRX(2, :, :))*O(4, :)')' - (squeeze(sRX(3, :, :))*O(3, :)')';
        G(2, :) = G(2, :) + (squeeze(sRX(3, :, :))*O(2, :)')' - (squeeze(sRX(1, :, :))*O(4, :)')';
        G(3, :) = G(3, :) + (squeeze(sRX(1, :, :))*O(3, :)')' - (squeeze(sRX(2, :, :))*O(2, :)')';
        
        G = G/n;

        cost = sum(log(diag(cHp))) - mean(log(s))*cdim;
    end


end

