function X = NRSfM_Consensus(D, flag_regu, regu_type, regu_order)
% function X = NRSfM_Consensus(D)
%
% Solve NRSfM by obtaining consensus from part reconstructions
%
% Inputs:
%     D: Input 2D trajectory data (k-1 x p x f) nooo
%     D: Input 2D trajectory data (k x p x f) but with 3 axis 0s
%     W: Weights - observed (true) or missing (false) (3 x p x n)
% 
%
% Outputs:
%     X: 3D reconstruction                            (k x p x f)

rotK_ratio = 1-1e-5;
tID_total = tic;

% preprocessing
% D(3, 1) = 0;
D = pout_trans(D);



%% 1) RANDOM SAMPLING
tic;
% select random groups

% parameters to change for user
nsamp = 10; %4; % no less 3
mgroup = 50; %12;
lambda = 0.1;

idx = select_idx(D(1:2, :, :), nsamp, mgroup, lambda); % [55, 358 grupos]
disp(['select_idx: ' num2str(toc)]);

ngroup = size(idx, 2); %set(gcf,'color','w') backgruond of plot to white
% spy(idx);title("IDX of 365 ngruops of walking dataser") ;

%% cluster plot

% flag_imgplot_list = [1, 20, 50, 70];
% 
% coger la funcion de cluster 3d plot y hacerlo en 2d plots.
% quiero ploterar todas las prates 2d divididas...
% plotear las partes reconst ?
% plotear antes de consensu -> 3d plot
% 
% 
% for j=1:size(flag_imgplot_list,2)
%     for i=1:ngroup % reconstruye grupo por grupo!!!
%         plot_NRSfM([], D(:, idx(:, i), :), ...
%             D(:, idx(:, i), :), ...
%             [], flag_imgplot_list(j) );
%         
%     end
%     cluster_3d_plots();    
% end
% cluster = D(:, idx(:, 1), :);
% for i=1:size(cluster,3)
%     scatter3(cluster(1, :, i), cluster(3, :, i), cluster(2, :, i), 'b.');
%     axis equal; title(dataname); drawnow;
% end
% 
% 




%% 2) WEAK RECONSTRUCTION

tic;
% solve for each group
tID = tic;
Xi = cell(1, ngroup); % x grupos, cada grupo hay como 10 tray??
% CHANGE VARIABLES:
%--------------------
max_ite = 500; %<---
order_L = regu_order; 
flag_soft = true; % if false-> hard, if true soft
if regu_type == "HARD"; flag_soft = false; end

%--------------------
for i=1:ngroup % reconstruye grupo por grupo!!!
    Xi{i} = reconstruct(D(:, idx(:, i), :), rotK_ratio, max_ite, order_L, flag_soft, flag_regu);

    % enviamos las 10 trayctorias/puntos para todos los frames.
    if toc(tID) > 1
        disp(['reconstruct ' num2str(i) ' / ' num2str(ngroup)]);
        tID = tic;
    end
end

disp(['reconstruct total: ' num2str(toc)]);


%% 3) REFLECTION CORRECTION

tic;
% reflection correction between groups
r = part_reflection(Xi, idx);
for i=1:ngroup
    if r(i) < 0
        Xi{i}(3, :, :) = -Xi{i}(3, :, :);
    end
end
disp(['part_reflection: ' num2str(toc)]);

%% 4) CONSENSUS
tic;
% combine weak reconstructions
X = combine(D, Xi, idx);
disp(['combine: ' num2str(toc)]);

toc(tID_total);
fprintf(1, 'Total time: %dm %ds\n', ...
     int16( toc(tID_total)/60), int16(mod( toc(tID_total), 60)) );

end

