load('./Data/feature_trajectories_L.mat', 'feature_trajectories')

load('./Data/back_L_3D.mat', 'track3D')

load('./Data/feature_trajectories_L.mat')

load('./Data/back_L_3D.mat')




 W2d=[];

 GT.shape3D=[];

 list=[];


%%


% indexs=[];

 for i=1:size(track3D,3)     
%     for j=1:size(track3D,2)
%         ag=find(track3D(:,j,i)==[0; 0; 0]);
%    if ~isempty(ag)
%        indexs=[indexs j];
%    end
%    
%     end
% 
% end
% 
% %total_track3D=track3D;
% track3D(:,indexs,:)=[];
% for i=1:150
%     
   W2d=[W2d; track3D(1,:,i); track3D(2,:,i)] ;
   GT.shape3D=[GT.shape3D; track3D(1,:,i); track3D(2,:,i); track3D(3,:,i)] ; %size 450 x245
% %%plot(feature_trajectories(2,:,i)./3,feature_trajectories(1,:,i)./3,'.r');
% % plot3(total_track3D(1,:,i),total_track3D(2,:,i),total_track3D(3,:,i),'.b');
% % hold on
% % plot3(track3D(1,:,i),track3D(2,:,i),track3D(3,:,i),'.r');
% % %hold off
% %       axis equal
% % %      pause(0.1);
% %        pause
%     
 end

 %% plot 2d back
Dd = GT(1:2,:,:);
frames_num = size(Dd,3);
for i=1:frames_num
    
    plot(Dd(1,:,i),Dd(2,:,i),'.r');
    axis equal
    pause(0.1);

end




 %%
seq = 'back_sparse';
X = W2d; % F per points
n_frames = size(X,1)/3;
x = X(1:2:end, :);
y = X(2:2:end, :);
% z = X(3:3:end, :);

X_new = zeros(2, size(X,2), n_frames);
for i = 1:n_frames
%    X_new(:,:,i) = [x(i,:);y(i,:);z(i,:)];
   % in the paper, the author change z as y
   X_new(:,:,i) = [x(i,:);y(i,:)]; 
end

X = X_new;
save(['./Data/' seq '_rearranged.mat'], 'X');
