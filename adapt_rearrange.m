function X_new = adapt_rearrange(X)
% maybe think other way to do it, more efficiently ?
n_frames = size(X,1)/3;
x = X(1:3:end, :);
y = X(2:3:end, :);
z = X(3:3:end, :);

X_new = zeros(3,size(X,2),n_frames);
for i = 1:n_frames
%    X_new(:,:,i) = [x(i,:);y(i,:);z(i,:)];
   % in the paper, the author change z as y
   X_new(:,:,i) = [x(i,:);z(i,:);y(i,:)];
end

end