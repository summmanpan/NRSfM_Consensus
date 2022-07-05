clear all

load('./Data/feature_trajectories_L.mat', 'feature_trajectories')

load('./Data/back_L_3D.mat', 'track3D')
%%
W2d=[];
GTshape3D=[];
list=[];
indexs=[];

for i=1:size(track3D,3)
    for j=1:size(track3D,2)
       ag=find(track3D(:,j,i)==[0; 0; 0]);
       if ~isempty(ag)
           indexs=[indexs j];
       end
    end
end

total_track3D=track3D;
track3D(:,indexs,:)=[];

%%
W2d = zeros(2,size(track3D,2)+size(total_track3D,2),150);
for i=1:150
    W2d(1, 1:size(track3D,2),i) = track3D(1,:,i);
    W2d(1, size(track3D,2)+1 : end, i) = total_track3D(1,:,i);
    W2d(2, 1:size(track3D,2),i) = track3D(2,:,i);
    W2d(2, size(track3D,2)+1 : end, i) = total_track3D(2,:,i);
   
%     W2d(:,1:size(track3D,2),:) = track3D(1,:,i);
%     W2d(:,1:size(total_track3D,2),:) = track3D(1,:,i);
%     xs = [track3D(1,:,i); total_track3D(1,:,i)];
%     ys = [track3D(2,:,i); total_track3D(2,:,i)];
%     W2d=[W2d; xs ; ys] ; % ahora no funciona pq nº de col es diff, hay q
%     pensar como concatenar tanto los puntos azules como los rojos
    
%     GTshape3D=[GTshape3D; track3D(1,:,i); track3D(2,:,i); track3D(3,:,i)];

    % plot 2D
%      plot(feature_trajectories(2,:,i)./3,feature_trajectories(1,:,i)./3,'.r');
    % plot 3D
%      plot3(total_track3D(1,:,i),total_track3D(2,:,i),total_track3D(3,:,i),'.b');
%      hold on
%      plot3(track3D(1,:,i),track3D(2,:,i),track3D(3,:,i),'.r');
%      hold off
%      axis equal
%      pause(0.1);
%      pause
%      drawnow limitrate;
 end

%% data arrangement:

n_frames = size(W2d,1)/2;
x = W2d(1:2:end, :);
y = W2d(2:2:end, :);

W_back = zeros(2, size(W2d,2), n_frames);
for i = 1:n_frames
   % in the paper, the author change z as y
   W_back(:,:,i) = [x(i,:);y(i,:)]; % put them last tensor position the index of the frame
end
%% plot

for i = 1:n_frames
    plot(W_back(1,:,i), W_back(2,:,i),'.k')
    
    axis equal
    pause(0.1);
    pause
%     drawnow limitrate;
end

%% save
D = W_back;
save('./Data/back_sparse_rearranged','D')


%% ACTRIZ

load './Data/point2D_ref_actriz.mat'
W=[];
for i=1:102
    W=[W; AllShapes(1:68,i)'; AllShapes(69:136,i)'];
end

%%
W2d=W;
n_frames = size(W2d,1)/2;
x = W2d(1:2:end, :);
y = W2d(2:2:end, :);

W_actriz = zeros(2, size(W2d,2), n_frames);
for i = 1:n_frames
   % in the paper, the author change z as y
   W_actriz(:,:,i) = [x(i,:);y(i,:)]; % put them last tensor position the index of the frame
end

%%
for i = 1:n_frames
    plot(W_actriz(1,:,i), -W_actriz(2,:,i),'.k')
    
    axis equal
%     pause(0.1);
%     pause
    drawnow limitrate;
end

%%
D = W_actriz;
save('./Data/actriz_rearranged','D')

