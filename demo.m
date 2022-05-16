% Demo program for NRSfM_Consensus.

% NRSfM_Consensus is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Loading GT 3D data
clear; close all; clc;

% GT should be (3 x p x f) dimensional, of the form;
% To run it, we only need the annotations, ie. GT(1:2,:,:). With following
% forms:
% GT(1:2,:,:) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% GT(:, :, k) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% Here, x_ki, y_ki, and z_ki are the 3D coord. of the ith point of the kth frame.
% W is the weight mask that indicates whether the
% elements are missing or not. (false = missing)


dataname = 'yoga'; % choice of data set 
% drink, pickup, stretch, yoga

% datatype = 'benchmark/';
% datatype = 'symthetic_camera_rotations/';

datapath = './Data/';
filename = '_rearranged.mat';

data = [datapath,dataname,filename];
load(data) % load the X as GT,or save it as X
GT = X;

%% Dense Data Set
% seq = {'nikos', 'back', 'heart'}; % mostñy imposible run in my computer...
ss  = 1;
dataname = seq{ss};
load(['./Data/dense/' seq{ss} '_rearranged.mat']);

%% Input data generation

% Experimental setting
noise = 0;            % noise level. paper use 10^-3
rmiss = 0;            % missing rate, value lower than 1

[k, p, nSample] = size(GT);
D = zeros(k, p, nSample);
temp = GT(1:2, :, :);
% add normal random number
weight_noise = noise*max(abs(reshape(bsxfun(@minus, temp, mean(temp, 2)), [], 1)));
D(1:2, :, :) = temp + weight_noise*randn(2, p, nSample);

% missing data
W = true(k, p, nSample);
W(3, :, :) = false;

ind = rand(p*nSample, 1) < rmiss;
D(1:2, ind) = 0; % set them to 0
W(1:2, ind) = false;

% Diary save for command window
resulpath = './Results_error/with_paper_rot/';
final_name_data = sprintf("%s%s",resulpath,'RESULTS_',dataname,'_120_2.txt');
% diary(final_name_data)
% d = datetime(now,'ConvertFrom','datenum'); % remove
disp(['**********************'+string(datetime)+'**********************'])
disp(['----------------------'+string(dataname)+'----------------------'])

% Consensus of Non-Rigid Reconstructions
X = NRSfM_Consensus(D,W);

%
% load Reconst_matlab_files\X_yoga.mat

% Evaluation
GT = bsxfun(@minus, GT, mean(GT, 2));
vind = sum((GT(3, :, :)-X(3, :, :)).^2) > sum((GT(3, :, :)+X(3, :, :)).^2);
% > q hay mas diff GT que X, es decir si GT es 10, y X 7, pues la diff si
% es mas mas grande que la suma de 10 mas 7 , quiero decir que
% si es 10-10=0 > 20-> no errror
% si error es mayor que 0, o no error, lo invertimos.
% si es 3 > ?? 1-10=-9 > 10 -> falso-> es neg -> invertimos
% aseguramos de que la diff de GT - X, es positivo 

X(3, :, vind) = -X(3, :, vind);

perf = sqrt(reshape(sum(sum((GT-X).^2)), 1, [])./reshape(sum(sum(GT.^2)), 1, []));
disp("--------------------------MEAN ERROR---------------------------")
disp(['mean error : ' num2str(mean(perf))]); % string(dataname)+':'+

% diary off

%% Plot 3D results

% for i=1:numel(perf)
%     scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');% LO TIENE REVES, Y PINTA Y -> AXIZ Z
%     hold on; scatter3(X(1, :, i), X(3, :, i), X(2, :, i), 'b.'); hold off;
%     axis equal; title(dataname); drawnow;
% end


% v_obj = VideoWriter(['./results/videos/' dataname '_video.avi']);
% plot_NRSfM(D, W, GT, X, plot_NRSfM(D, W, GT, X););

% list = [];
plot_NRSfM(D, W, list, GT, X);



%% Save VARIABLES, the reconstruct X
% X_path = ['./Reconst_matlab_files/X_'+string(dataname)+'.mat'];
% save(X_path,"X")
% saveas(Figure 1,dataname)

% clusters_3d_plots(Xi)
