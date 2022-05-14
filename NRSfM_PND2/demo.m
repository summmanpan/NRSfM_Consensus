% Demo program for NRSfM_PND2.
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



% **IMPORTANT: Run 'download_and_rearrange_data.m' at the first time to download data sets.
% If you are going to use another data set for this demo, the data should be rearranged as
%
% X(:, :, k) = [x_k1 x_k2 ... x_kp; y_k1 y_k2 ... y_kp; z_k1 z_k2 ... z_kp]
%
% Here, x_ki and y_ki are the observed 2D coordinates of the kth frame, and z_ki are the
% corresponding depth coordinates. Then, the observation data D will be generated as
% D(:, :, k) = [X(1:2, :, k); zeros(1, p)]. W is the weight mask that indicates whether the
% elements are missing or not. (false = missing)


clear; close all; clc;


% Experimental setting
seq = 'walking';    % choice of data set
noise = 0.0001;          % noise level
rmiss = 0.001;          % missing rate


% Loading
load(['Data/' seq '_rearranged']);
[k, p, nSample] = size(X);

% Input data generation
D = zeros(k, p, nSample);
temp = X(1:2, :, :);
D(1:2, :, :) = temp + noise*max(abs(reshape(bsxfun(@minus, temp, mean(temp, 2)), [], 1)))*randn(2, p, nSample);

W = true(k, p, nSample);
W(3, :, :) = false;

ind = rand(p*nSample, 1) < rmiss;
D(1:2, ind) = 0;
W(1:2, ind) = false;

tic;

% NRSfM: EM-PND (Procrustean Normal Distribution)
[rX, s, R, M, C, var_D] = NRSfM_PND2(D, W);

toc;

% Evaluation
gX = bsxfun(@minus, X, mean(X, 2));
vind = sum((gX(3, :, :)-rX(3, :, :)).^2) > sum((gX(3, :, :)+rX(3, :, :)).^2);
rX(3, :, vind) = -rX(3, :, vind);

perf = sqrt(reshape(sum(sum((gX-rX).^2)), 1, [])./reshape(sum(sum(gX.^2)), 1, []));
disp(['mean error : ' num2str(mean(perf))]);

% Plot results
plot_NRSfM(D, W, gX, rX);

