function X = adapt_rearrange(X, seq)
% P3_gt : GT 3D information in matrix
% seq: string name of the seq
% trD: flag of it has Z axis

n_frames = size(X,1)/3;
x = X(1:3:end, :);
y = X(2:3:end, :);
z = X(3:3:end, :);

X_new = zeros(3,size(X,2),n_frames);
for i = 1:n_frames
%    X_new(:,:,i) = [x(i,:);y(i,:);z(i,:)];
   % in the paper, the author change z as y
   X_new(:,:,i) = [x(i,:);z(i,:);y(i,:)]; % put them last tensor position the index of the frame
end

X = X_new;
save(['./Data/' seq '_rearranged.mat'], 'X');

% rearanged it more efficiently:
%  Input: 
%  P3_gt (3f x p)
% 
%  Output:
%  X     (3 x p x f)
% axis=0;
% if triD == 1 % has axis 3
%     axis = 3;
% end
% [~,p] = size(P3_gt);
% X = reshape(reshape(P3_gt , [] , axis*p)', axis, p, []);

end