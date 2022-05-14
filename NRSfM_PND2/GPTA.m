function [Y, s, R, M] = GPTA(X)
% function [Y, s, R, M] = GPTA(X)
%
% Modified generalized Procrustes analysis
%
% Inputs:
%     X: Input shapes   (3 x p x n)
%
% Outputs:
%     Y: Aligned shapes (3 x p x n)
%     s: Scales         (1 x n)
%     R: Rotations      (3 x 3 x n)
%     M: Mean shape     (3 x p)
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

m = size(X, 3);

Y = bsxfun(@minus, X, mean(X, 2));
s = 1./reshape(sqrt(sum(sum(Y.^2))), 1, m);
Y = bsxfun(@times, Y, reshape(s, 1, 1, m));
M = mean(Y, 3);
M = M/norm(M(:));

count = 0;
pM = zeros(size(M));
R = zeros(3, 3, m);
R(1, 1, :) = 1;
R(2, 2, :) = 1;
R(3, 3, :) = 1;
while mse(M-pM) > eps
    pM = M;
    
    for i=1:m
        [U, S, V] = svd(Y(:, :, i)*M');
        ts = 1/trace(S);
        tR = V*U';
        s(i) = ts*s(i);
        R(:, :, i) = tR*R(:, :, i);
        Y(:, :, i) = ts*tR*Y(:, :, i);
    end
    
    M = mean(Y, 3);
    M = M/norm(M(:));

    count = count + 1;
%     disp([num2str(count) ' : ' num2str(mse(M-pM)) ' / ' num2str(mse(bsxfun(@minus, bsxfun(@rdivide, reshape(Y, 3*n, m), M(:)'*reshape(Y, 3*n, m)), M(:))))]);
end

[U, ~, ~] = svd(M, 'econ');
M = U'*M;
Y(:, :) = U'*Y(:, :);
R(:, :) = U'*R(:, :);

end

