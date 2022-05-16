function X_new = adapt_rearrange(P3_gt, triD)

% n_frames = size(X,1)/3;
% x = X(1:3:end, :);
% y = X(2:3:end, :);
% z = X(3:3:end, :);
% 
% X_new = zeros(3,size(X,2),n_frames);
% for i = 1:n_frames
% %    X_new(:,:,i) = [x(i,:);y(i,:);z(i,:)];
%    % in the paper, the author change z as y
%    X_new(:,:,i) = [x(i,:);z(i,:);y(i,:)]; % put them last tensor position the index of the frame
% end

% rearanged it more efficiently:
%  Input: 
%  P3_gt (3f x p)
% 
%  Output:
%  X     (3 x p x f)
axis=0;
if triD == 1
    axis = 3;
end
[~,p] = size(P3_gt);
X = reshape(reshape(P3_gt , [] , axis*p)', axis, p, []);
save(['Data_rearranged/' seq '_rearranged.mat'], 'X');
disp([seq ' done.']);
end