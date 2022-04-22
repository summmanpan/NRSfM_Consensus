% Demo program for NRSfM_Consensus.

%
% NRSfM_Consensus is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

clear; close all; clc;

%*************************************************************************************
% TODO: You have to load a ground truth 3D data (named "GT") here, e.g.;

datapath = './data_set/';
datatype = 'benchmark/';
dataname = 'walking'; %'walking';

% datatype = 'symthetic_camera_rotations/';
% dataname = 'drink'; % drink, pickup, stretch, yoga
filename = '_rearranged.mat';
data = [datapath,datatype,dataname, filename];

% data = [datapath,dataname, filename];
load(data)

% GT should be (3 x p x f) dimensional, of the form;
%  La matriz de entrada es 3 X puntos X imágenes.
%  Para ejecutarlo solo necesitas la W, las anotaciones en la imagen, 
% es decir, GT(1:2,:,:).
GT = X;
% GT(1:2,:,:) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% GT(:, :, k) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];

% Here, x_ki, y_ki, and z_ki are the 3D coordinates of the ith point of the kth frame.
%*************************************************************************************
%% add rotation
ang = 60;
GT_rot = addRotation(ang,GT,dataname); 
% intentar de crear que pueda add y remove rotation entre un rg de frames
% que quiera el usuario...
% creo q deberia funcionar con el inverse matrix, pero no ! mirar por que!
plot_2D(GT_rot,dataname)
%%
% plot_2D(GT,dataname)
% Input data generation
D = GT(1:2, :, :);

% Diary save for command window
resulpath = './Results_error/with_paper_rot/';
final_name_data = sprintf("%s%s",resulpath,'RESULTS_',dataname,'_120_2.txt');
% diary(final_name_data)
d = datetime(now,'ConvertFrom','datenum');
disp(['**********************'+string(datetime)+'**********************'])
disp(['----------------------'+string(dataname)+'----------------------'])

% Consensus of Non-Rigid Reconstructions
X = NRSfM_Consensus(D);

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

for i=1:numel(perf)
    scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');% OJOO Q EL LO TIENE REVES, Y PINTA Y -> AXIZ Z
    hold on; scatter3(X(1, :, i), X(3, :, i), X(2, :, i), 'b.'); hold off;
    axis equal; title(dataname); drawnow;
end

% for i=1:numel(perf)
%     clf;
%     scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
%     hold on; scatter3(X(1, :, i), X(3, :, i), X(2, :, i), 'b.'); hold off;
%     axis equal;
%     hold on;
%     title(dataname);
%     getframe;
%     pause;
% end

%% Save the reconstruct X
% X_path = ['./Reconst_matlab_files/X_'+string(dataname)+'.mat'];
% save(X_path,"X")
% saveas(Figure 1,dataname)

%% 3D plots Xi reconstruct parts CLUSTERS
hFig = figure();
axh = axes('Parent', hFig);
hold(axh, 'all');
grid(axh, 'on');
frames_number = size(Xi{1}(:,:,:),3);
for i=1: frames_number
    number_end =15;
    for s = 1:2
        if s == 1
            subplot(1,3,1)
        else
            subplot(1,3,2)
        end
        for j = 1:number_end %
            r_part = Xi{j}(:,:,:); 
            scatter3(r_part(1, :, i), r_part(3, :, i), r_part(2, :, i),'filled');
            hold on
            if j==number_end % llega al final
                hold off
            end
        end
       
        
        if s == 1
           %subplot(1,3,3)
%            hold on
%            scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
           view(3);
%            axis equal; drawnow;
           %view([AZ,EL])
           hold off
        else
            view(2); 
        end
        axis equal; drawnow;
    end
%     subplot(1,3,3)
%     scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
%     view(3);
%     axis equal; drawnow;  
end

% legend(axh, [h1,h2], {'Alpha', 'Beta'});
xlabel('X')
ylabel('Y')
zlabel('Z')

%%

R*S(3i-2:3i,:); % matrix de size 3xp
S_new = [];

% for frames
% Las 1 y 2 partes -> hacerlos que sean multiplicados con R=I? -> pensar y
% hay otra forma, o simplemente hacerlos con las partes que tocan
S_new = [S_new; R*S(3i-2:3i,:)]; % esto los contatena ; size 3fxp
% los posteriores tamb

% hacer una comparacion del plot , con rotada y sin.
% puedo usar la variable lista para hacer las lineas y tener una mejor
% visualizacion!

% luego adaptar el data al de Lee!


