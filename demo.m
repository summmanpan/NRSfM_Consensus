% Demo program for NRSfM_Consensus.
%
% Ref: Minsik Lee, Jungchan Cho, and Songhwai Oh,
% "Consensus of Non-Rigid Reconstructions,"
% CVPR 2016, Las Vegas, Nevada, June 26-July 1, 2016.
%
% Author: Minsik Lee (mlee.paper@gmail.com)
% Last update: 2016-08-08
% License: GPLv3

%
% Copyright (C) 2016 Minsik Lee
% This file is part of NRSfM_Consensus.
%
% NRSfM_Consensus is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% NRSfM_Consensus is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with NRSfM_Consensus. If not, see <http://www.gnu.org/licenses/>.


clear; close all; clc;


%*************************************************************************************
% TODO: You have to load a ground truth 3D data (named "GT") here, e.g.;
% 
% load dance; % This file should contain GT
load ./data/yoga_rearranged.mat
%
% 
% GT should be (3 x p x f) dimensional, of the form;
% p 
% el dataset tendrá x, y, z, dimension 3D
%
%  La matriz de entrada es 3 X puntos X imágenes.
%  Para ejecutarlo solo necesitas la W, las anotaciones en la imagen, 
% es decir, GT(1:2,:,:).
GT = X;
% GT(1:2,:,:) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% GT(:, :, k) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
    
%
% Here, x_ki, y_ki, and z_ki are the 3D coordinates of the ith point of the kth frame.
%*************************************************************************************

%%

% Input data generation
D = GT(1:2, :, :);

% Consensus of Non-Rigid Reconstructions
X = NRSfM_Consensus(D);

% Evaluation
GT = bsxfun(@minus, GT, mean(GT, 2));
vind = sum((GT(3, :, :)-X(3, :, :)).^2) > sum((GT(3, :, :)+X(3, :, :)).^2);
X(3, :, vind) = -X(3, :, vind);

perf = sqrt(reshape(sum(sum((GT-X).^2)), 1, [])./reshape(sum(sum(GT.^2)), 1, []));
disp(['mean error : ' num2str(mean(perf))]);

% Plot results
for i=1:numel(perf)
    scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
    hold on; scatter3(X(1, :, i), X(3, :, i), X(2, :, i), 'b.'); hold off;
    axis equal; drawnow;
end

