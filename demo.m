% Demo program for NRSfM_Consensus.

% GT should be (3 x p x f) dimensional, of the form;
% To run it, we only need the annotations, ie. GT(1:2,:,:).
% GT(1:2,:,:) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% GT(:, :, k) = [ x_k1 x_k2 x_kp; y_k1 y_k2 y_kp; z_k1 z_k2 z_kp];
% Here, x_ki, y_ki, and z_ki are the 3D coord. of the ith point of the kth frame.

% W is the weight mask that indicates whether the
% elements are missing or not. (false = missing)

%% Dense Data Set
% seq = {'nikos', 'back', 'heart'}; % mostly imposible to run in my computer...
% ss  = 2;
% dataname = seq{ss};
% load(['./Data/dense/' seq{ss} '_rearranged.mat']);
% 
% if size(GT,1)==2
%     D=GT;
%     D(3,1)=0;
% end


%% Automatizar el testing

% clear; close all; clc;

with_noise = {'Noise','NoNoise'};
with_regu= {'Regu','NoRegu'};

regu_type = {'SOFT','HARD'}; 
regu_order = [1,2,4]; 
with_rot = {'Rot','noRot'};
rot_list = {'120'}; % '','60','90','120' " % save for each dataset rotation individually<-------

dataname = {'drink','pickup', 'stretch', 'yoga'};  %drink   
% dataname = {'drink'};  %drink   ESTE COMO TARDA MUCHO, LO PONDRE A PARTE!!
%     dataname = {'drink','pickup', 'stretch', 'yoga'}; % drink -> ESTE APARTE PERO LUEGO PONER!!!
%     dataname = {'walking','dance','jaws','face'};

%-------------------- CHANGE PARAMETRES------------------------
flag_noise = 2 ; %<--
flag_regu = 2 ;% 1- regu, 2-no Regu
% only when with regu:
flag_regu_type = 2 ;% 1-soft or 2-hard

flag_rot = 1 ; % 1-rot 2-no rot
if flag_rot ==2; rot_list = {''}; end
noise = 0.001;
odr = 3;


% for type_r = 1:size(regu_type,2)
    
%     for odr = 1 :size(regu_order,2)
        
        X_cell = cell(size(rot_list,2) ,size(dataname,2) );
        err_cell = cell(size(rot_list,2) ,size(dataname,2) );
    
        for i=1:size(rot_list,2) % rot number
        
            for j= 1 : size(dataname,2) 
                
                % Load the data
                path_data = ['./Data/' dataname{j} '_rearranged.mat'];
                if rot_flag==1
                    path_data = ['./Data/rot_' rot_list{i} '/' dataname{j} '_' rot_list{i} '_rearranged.mat'];
                end
                load(path_data, 'X');
                
                % Diary save for command window
                resulpath = './Results_error/text_error/';
                final_name_data = [resulpath 'ERROR_' dataname{j} '_' rot_list{i} '.txt'];
                
%                 diary(final_name_data)
                disp(['***'+string(datetime)+'***'])
                disp(['---'+string(dataname{j})+'---'])

%                 [X_cell{i,j}, err_cell{i,j}] = get_principal_function(X,dataname{j}, ...
%                                             flag_regu, regu_type{flag_regu_type}, ...
%                                             regu_order(odr), noise); 
                
  
               
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
%                  diary off
            end
    
            principal_path= './Results_error/';
            noise_path = sprintf('%s/', with_noise{flag_noise});
            if flag_regu == 2  % no regu
                regu_path = sprintf('%s/data_%s_%s_%s', ...
                                    with_regu{flag_regu}, ...
                                    with_regu{flag_regu}, with_rot{flag_rot}, ...
                                    rot_list{i});
            else
                regu_path = sprintf('%s/%s/L%d/data_%s_%s_L%d_%s_%s,', ...
                                    with_regu{flag_regu},regu_type{flag_regu_type}, ...
                                    regu_order(odr), ... % data_final_path
                                    with_regu{flag_regu},regu_type{flag_regu_type}, ...
                                    regu_order(odr),with_rot{flag_rot}, ...
                                    rot_list{i});
            end
            final_save_path = [principal_path noise_path regu_path];
%             save(final_save_path,'X_cell','err_cell','noise')
            disp(final_save_path)
           
        end
%     end
% end



%% Plot 3D results

% for i=1:numel(perf)
%     scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');% LO TIENE REVES, Y PINTA Y -> AXIZ Z
%     hold on; scatter3(X(1, :, i), X(3, :, i), X(2, :, i), 'b.'); hold off;
%     axis equal; title(dataname); drawnow;
% end

% 
% v_obj = VideoWriter(['./results/videos/' dataname '_video.avi']);
% % plot_NRSfM(D, W, GT, X, plot_NRSfM(D, W, GT, X););
% 
% list = [];
% plot_NRSfM(D, W, list, GT, X, v_obj);



%% Save VARIABLES, the reconstruct X
% X_path = ['./Reconst_matlab_files/X_'+string(dataname)+'.mat'];
% save(X_path,"X")
% saveas(Figure 1,dataname)

% clusters_3d_plots(Xi)


% dataname = {'dinosaur_real','face_mocap','face_real','face','FRGC'}; %ogre_synthetic
%             if type_regu~=3
%                 path_save = sprintf('./Results_error/%s/L%d/data_%s_L%d_%s_rot_%s', ...
%                                             regu_type{type_regu}, regu_order(odr), ...
%                                             regu_type{type_regu}, regu_order(odr), ...
%                                              datatype{datatype_idx}, rot_list{i});
%                 save(path_save,'X_cell','err_cell','noise')
%             end
            
%             path_save = sprintf('./Results_error/%s/data_%s_%s_rot_%s_Motion', ...
%                         regu_type{type_regu},regu_type{type_regu}, ...
%                          datatype{datatype_idx},rot_list{i});
%             save(path_save,'X_cell','err_cell','noise')

                % save for each dataset rotation individually
                
%                 if dataname{1} == "drink"
%                     path_save = sprintf('./Results_error/%s/L%d/data_%s_L%d_%s_rot_%s', ...
%                                         regu_type{type_regu}, regu_order(odr), ...
%                                         regu_type{type_regu}, regu_order(odr), ...
%                                          dataname{1}, rot_list{i}  );
%                 end
%                 save(path_save,'X_cell','err_cell')


% datatype = {'originalRotation','noRotation', 'benchmark'};

% path_save = sprintf('./Results_error/datasets_rot_%s_withReguHardOrder_%d', rot_list{1}, regu_order(1));
% path_save = sprintf('./Results_error/datasets_rot_%s_withReguSoftOrder_%d', rot_list{1}, regu_order(1));
% save(path_save,'X_cell','err_cell')

