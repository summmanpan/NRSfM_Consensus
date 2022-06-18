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
for i=1:150
    W2d=[W2d; track3D(1,:,i); track3D(2,:,i)] ;
    GTshape3D=[GTshape3D; track3D(1,:,i); track3D(2,:,i); track3D(3,:,i)] ;
    % plot 2D
%     plot(feature_trajectories(2,:,i)./3,feature_trajectories(1,:,i)./3,'.r');
    % plot 3D
%     plot3(total_track3D(1,:,i),total_track3D(2,:,i),total_track3D(3,:,i),'.b');
%     hold on
%     plot3(track3D(1,:,i),track3D(2,:,i),track3D(3,:,i),'.r');
%     hold off
%     axis equal
%     pause(0.1);
%     pause

 end



% elimina una columna de %, y en teoria utilizando esto puedes tener la W2d de entrada

