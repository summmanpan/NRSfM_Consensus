function X = NRSfM_Consensus(D)
% function X = NRSfM_Consensus(D)
%
% Solve NRSfM by obtaining consensus from part reconstructions
%
% Inputs:
%     D: Input 2D trajectory data                     (k-1 x p x f)
%
% Outputs:
%     X: 3D reconstruction                            (k x p x f)

rotK_ratio = 1-1e-5;


tID_total = tic;

% preprocessing
D(3, 1) = 0;
D = pout_trans(D);

tic;
% select random groups
nsamp = 10;
mgroup = 50;
lambda = 0.1;
idx = select_idx(D(1:2, :, :), nsamp, mgroup, lambda);
disp(['select_idx: ' num2str(toc)]); % why we print the toc ? 

ngroup = size(idx, 2);

tic;
% solve for each group
tID = tic;
Xi = cell(1, ngroup); % x grupos, cada grupo hay com0 10 tray??
for i=1:ngroup
    Xi{i} = reconstruct(D(:, idx(:, i), :), rotK_ratio);
    if toc(tID) > 1
        disp(['reconstruct ' num2str(i) ' / ' num2str(ngroup)]);
        tID = tic;
    end
end
disp(['reconstruct total: ' num2str(toc)]);

tic;
% reflection correction between groups
r = part_reflection(Xi, idx);
for i=1:ngroup
    if r(i) < 0
        Xi{i}(3, :, :) = -Xi{i}(3, :, :);
    end
end
disp(['part_reflection: ' num2str(toc)]);

tic;
% combine weak reconstructions
X = combine(D, Xi, idx);
disp(['combine: ' num2str(toc)]);

toc(tID_total);

end
