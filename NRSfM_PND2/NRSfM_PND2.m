function [X, s, R, M, C, var_D] = NRSfM_PND2(D, W)
% function [X, s, R, M, C, var_D] = NRSfM_PND2(D, W)
%
% Solve NRSfM problem based on PND2
%
% Inputs:
%     D: Observations                                 (3 x p x n)
%     W: Weights - observed (true) or missing (false) (3 x p x n)
%
% Outputs:
%     X: Reconstructed shapes (3 x p x n)
%     s: Scales               (1 x n)
%     R: Rotations            (3 x 3 x n)
%     M: Mean shape           (3 x p)
%     C: Shape covariance     (3p x 3p)
%     var_D: Noise variance
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


% initialization
[D, nD, s, R, X, M] = initialize(D, W); % <- HABRA QUE MIRAR LAS W 

% EM-PND2
[X, s, R, M, C, var_D] = EM_PND2(s, R, X, M, D, W);%<- Y AQUI

% Restore the scale
X = X*nD;
s = s/nD;
var_D = var_D*nD^2;

end

