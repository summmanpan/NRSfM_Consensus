
clc ;
clear all;
% load the data
seq = {'nikos', 'back', 'heart'};
ss  =3;
load(['./Data/dense/' seq{ss} '_dense.mat']);

T=size(W2D,1)/2; % number of frames
W2d=zeros(size(W2D));

for i=1:T

    W2d(2*i-1,:)=W2D(i,:) ; 
    W2d(2*i,:)=W2D(i+T,:); 

end
%%
% for i = 1:T
%     plot(W2d(2*i-1,:), W2d(2*i,:),'.k','MarkerSize',1)
%     drawnow
% end

n_frames = size(W2d,1)/2;
x = W2d(1:2:end, :);
y = W2d(2:2:end, :);

GT = zeros(2,size(W2d,2),n_frames);
for i = 1:n_frames
   % in the paper, the author change z as y
   GT(:,:,i) = [x(i,:);y(i,:)]; % put them last tensor position the index of the frame
end

for i = 1:T
    plot(GT(1,:,i), GT(2,:,i),'.k','MarkerSize',1)
    drawnow
end

save(['./Data/dense/' seq{ss} '_rearranged.mat'], 'GT');

