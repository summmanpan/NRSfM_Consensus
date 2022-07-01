% Demo program for NRSfM_Consensus.

% GT should be (3 x p x f) dimensional, of the form;
% To run it, we only need the annotations, ie. GT(1:2,:,:).
% GT(1:2,:,:) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% GT(:, :, k) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% Here, x_ki, y_ki, and z_ki are the 3D coord. of the ith point of the kth frame.

% W is the weight mask that indicates whether the
% elements are missing or not. (false = missing)
%% Sparse back 

dataname = 'back_sparse'; % mostly imposible to run in my computer..., 'back'

load(['./Data/' dataname '_rearranged.mat']);

if size(D,1)==2
    D(3,1)=0;
end
%%
flag_regu = 1; % 1 si regu, 2 no regu, 
[X, ~] = get_principal_function(D, dataname, ...
                    flag_regu, 'HARD', ...
                    3, 0); 

%% Dense Data Set
% dataname = {'back'}; % mostly imposible to run in my computer..., 'back'
% flag_regu=1;
% for ss = 1:1 %size(seq,2)
%     load(['./Data/dense/' dataname{ss} '_rearranged.mat']);
%     
%     if size(GT,1)==2
%         D=GT(:,1:40:end,:);
%         D(3,1)=0;
%     end
%     [X, ~] = get_principal_function(D, dataname{ss}, ...
%                         flag_regu, 'SOFT', ...
%                         1, 0); 
%     % mñana compiar lo de abajo y correr para guardarlo
% %     path_save=sprintf("./results/data_%s_Regu_SOFT_1_NoNoise",dataname{ss});
% %     save(path_save)
% end

%% Automatizar el testing

% clear; close all; clc;

with_noise = {'Noise','NoNoise'};
with_regu= {'Regu','NoRegu'};

regu_type = {'SOFT','HARD'}; 
regu_order = [1,2,4]; 

with_rot = {'Rot','noRot'};
rot_list = {'90'}; % '','60','90','120' " % save for each dataset rotation individually<-------

% dataname = {'pickup', 'stretch', 'yoga'};  %drink   
% dataname = {'yoga'};  %drink   ESTE COMO TARDA MUCHO, LO PONDRE A PARTE!!
%     dataname = {'drink','pickup', 'stretch', 'yoga'}; % drink -> ESTE APARTE PERO LUEGO PONER!!!
%     dataname = {'walking','dance','jaws','face'};
dataname = {'dance'};

%-------------------- CHANGE PARAMETRES------------------------
flag_noise = 2 ; % 1- noise, 2- no noise
flag_regu = 1 ; % 1- regu, 2-no Regu
flag_regu_type = 2;% 1-soft or 2-hard % only when with regu:
flag_rot = 2 ; % 1-rot 2-no rot
 
if flag_rot == 2; rot_list = {''}; end
noise = 0.001;
if flag_noise == 2; noise=0; end
% odr = 3;

% for flag_regu_type = 1:size(regu_type,2)
    
    for i=1:size(rot_list,2) % rot number
        
        for odr = 3 : 3%size(regu_order,2) 
            X_cell = cell(size(1,2) ,size(dataname,2) );
            err_cell = cell(size(1,2) ,size(dataname,2) );
    
            for j= 1 : size(dataname,2) 
%                 dataname, flag_rot,rot_list, flag_list
                [X,~] = get_load_data(dataname{j}, flag_rot,rot_list{i}, 0);
                
                
                % Diary save for command window
                resulpath = './Results_error/text_error/';
                final_name_data = [resulpath 'ERROR_' dataname{j} '_' rot_list{i} '.txt'];
                
                diary(final_name_data)
                disp(['***'+string(datetime)+'***'])
                disp(['---'+string(dataname{j})+'---'])

                [X_cell{1,j}, err_cell{1,j}] = get_principal_function(X,dataname{j}, ...
                                            flag_regu, regu_type{flag_regu_type}, ...
                                            regu_order(odr), noise); 
                
  
               
                if flag_noise ==1
                    disp(['***WITH:***''Noise rate: ' num2str(noise) '***'])
                end
                if flag_regu == 2 % NO regu
                    
                    str = sprintf('%s_%s_%s', ...
                                    with_regu{flag_regu}, with_rot{flag_rot}, ...
                                    rot_list{i});
                    
                else
                    str = sprintf('%s_%s_L%d_%s_%s', ...
                                    with_regu{flag_regu},regu_type{flag_regu_type}, ...
                                    regu_order(odr),with_rot{flag_rot}, ...
                                    rot_list{i});
                end
                disp(str)
                 % set diary off
                diary off
            end
    
            principal_path= './Results_error/';
            noise_path = sprintf('%s/', with_noise{flag_noise});
            if flag_regu == 2  % no regu
                regu_path = sprintf('%s/data_%s_%s_%s', ...
                                    with_regu{flag_regu}, ...
                                    with_regu{flag_regu}, with_rot{flag_rot}, ...
                                    rot_list{i});
            else
                regu_path = sprintf('%s/%s/L%d/data_%s_%s_L%d_%s_%s', ...
                                    with_regu{flag_regu},regu_type{flag_regu_type}, ...
                                    regu_order(odr), ... % data_final_path
                                    with_regu{flag_regu},regu_type{flag_regu_type}, ...
                                    regu_order(odr),with_rot{flag_rot}, ...
                                    rot_list{i});
            end

            final_save_path = [principal_path noise_path regu_path];
            if dataname{j}=="drink"
                final_save_path = [final_save_path 'drink'];
            elseif dataname{j}=="dinosaur_real"
                final_save_path = [final_save_path 'dinosaur'];
            end

%             save(final_save_path,'X_cell','err_cell','noise')
            disp(final_save_path)
           
        end
    end
% end



%% Plot 3D results

for i=1:size(GT,3)
    scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');% LO TIENE REVES, Y PINTA Y -> AXIZ Z
    hold on; 
    scatter3(rX(1, :, i), rX(3, :, i), rX(2, :, i), 'b.');
    hold off;
    axis equal; title(dataname); drawnow;
end


%% GENERATE VIDEO ETC

% GT
dataname = 'dance';
flag_rot=0;
rot_list ='90';
flag_list=0;
[GT,list] = get_load_data(dataname, flag_rot, rot_list, flag_list);

%% get reconstruct data from folder


% with_noise = {'Noise','NoNoise'};
% with_regu= {'Regu','NoRegu'};

% regu_type = {'SOFT','HARD'}; 
% regu_order = [1,2,4]; 
% with_rot = {'Rot','noRot'};
% rot_list = {'120'}; % '','60','90','120' " % save for each dataset rotation individually<-------
% 
% dataname = {'drink','pickup', 'stretch', 'yoga'};  
% dataname = {'walking','dance','jaws','face'};

%-------------------- CHANGE PARAMETRES------------------------
% flag_noise = 1 ; % 1- noise, 2- no noise
% flag_regu = 1 ; % 1- regu, 2-no Regu
% 
% flag_regu_type = 1 ;% 1-soft or 2-hard % only when with regu:
% flag_rot = 1 ; % 1-rot 2-no rot
%  
rX = X_cell{1,2};
%%

list = [];
% rX = X_waliking;
% extras = [rot_list '__SOFT_L4'];
% video saver:
v_obj = VideoWriter(['./results/videos/' dataname '_back_C_H2_video.avi']);
plot_NRSfM(list, GT, rX, v_obj,0);

%% Export figure frame by frame
list = [];
% rX = X;
flag_imgplot_list = [25,45,65,90];%[230, 110, 150,188, 25, 50,75, 100,150,200,250,300];
for i=1:size(flag_imgplot_list,2)
    plot_NRSfM(list, GT, rX, [], flag_imgplot_list(i) );
end

% EN LA PRESENTATTION YA PLOTEO LOS CLUSTERS TAMBIEN!

%% Save VARIABLES, the reconstruct X
% X_path = ['./Reconst_matlab_files/X_'+string(dataname)+'.mat'];
% save(X_path,"X")
% saveas(Figure 1,dataname)
% clusters_3d_plots(Xi)

%%%%%%%%%%%%%%%%%%%%%5
% 
% function clickonBody(){
%     console.log("Clicked on body"); 
%     document.body.click()
% }
% setInterval(clickonBody,600000)
% Load the data
                
%% UTILS

function [X,list] = get_load_data(dataname, flag_rot,rot_list, flag_list)

    path_data = ['./Data/' dataname '_rearranged.mat'];
    if flag_rot==1
        path_data = ['./Data/rot_' rot_list '/' dataname '_' rot_list '_rearranged.mat'];
    end
    list=[];
    if flag_list == 1
        path_list = ['./Data/list_points/list_' dataname '.mat'];
        load(path_list,'list');
    end
    load(path_data, 'X');

end