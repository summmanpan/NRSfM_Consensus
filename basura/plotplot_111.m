n_frames = size(W,1)/2;
x = W(1:2:end, :);
y = W(2:2:end, :);

X_new = zeros(2,size(W,2),n_frames);
for i = 1:n_frames
%    X_new(:,:,i) = [x(i,:);y(i,:);z(i,:)];
   % in the paper, the author change z as y
   X_new(:,:,i) = [x(i,:);y(i,:)]; % put them last tensor position the index of the frame
end

%%
% dataname = {'heart'}; % mostly imposible to run in my computer..., 'back'
dataname="face_r";
flag_regu=1;
GT = X_new;
for ss = 1:1 %size(seq,2)
%     load(['./Data/dense/' dataname{ss} '_rearranged.mat']);
    
    if size(GT,1)==2
        D=GT;
        D(3,1)=0;
    end
    [X, ~] = get_principal_function(D, dataname, ...
                        flag_regu, 'SOFT', ...
                        1, 0); 
    % mñana compiar lo de abajo y correr para guardarlo
%     path_save=sprintf("./results/data_%s_Regu_SOFT_1_NoNoise",dataname{ss});
%     save(path_save)
end