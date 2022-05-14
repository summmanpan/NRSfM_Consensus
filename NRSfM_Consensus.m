function X = NRSfM_Consensus(D,W)
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

%% 1) RANSOM SAMPLING
tic;
% select random groups
nsamp = 10; %must be less than 55
mgroup = 50;
lambda = 0.1;
idx = select_idx(D(1:2, :, :), nsamp, mgroup, lambda); %[55, 358 grupos]
disp(['select_idx: ' num2str(toc)]); % why we print the toc ? 

ngroup = size(idx, 2); %set(gcf,'color','w') backgruond of plot to white
% spy(idx);title("IDX of 365 ngruops of walking dataser") ;

% SI AL FINAL AÑADO EL MISSING DATA CON MATRIX COMPLETION
% TENGO QUE HACER W_index?? CREO Q SERA MUY COMPLICADO:(

%% 2) WEAK RECONSTRUCTION

tic;
% solve for each group
tID = tic;
Xi = cell(1, ngroup); % x grupos, cada grupo hay com0 10 tray??
for i=1:ngroup % reconstruye grupo por grupo!!!
    Xi{i} = reconstruct(D(:, idx(:, i), :), rotK_ratio);
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
