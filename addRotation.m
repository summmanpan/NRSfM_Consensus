function S_new = addRotation(ang,S,dataname,list)
% according to the paper:
% starting from the 1/4 position of the entire frames
% ending at th2 half position of the frame
% the camera rorations were made around the yaxis

% syms R_y(ang)s
% R_y(ang) = [[cos(ang) 0 sin(ang)];[0 1 0];[-sin(ang) 0 cos(ang)]];

total_frame = size(S,1)/3;
start_pos = floor(total_frame/4);
end_pos = floor(total_frame/2);
ang_per_frame = ang/(end_pos-start_pos); % solo el cuarto parte
rand_y = 0;
rad = deg2rad( ang_per_frame );
R_y=0;

% S_new=[];
S_new = S(1:start_pos*3,:);
for i=start_pos+1:end_pos % per frame
    rand_y = rand_y + rad;
%     ang_y = ang_y + ang_per_frame;
    % Ry = R_y(ang_per_frame);
%     R_y = [[cos(rand_y) 0 sin(rand_y)];[0 1 0];[-sin(rand_y) 0 cos(rand_y)]];
%     matrix_in = [-1 0 0; 0 1 0]
% R_y = [[cos(rand_y) -sin(rand_y) 0];[sin(rand_y) cos(rand_y) 0];[0 0 1]]; % minus axis Z
    % Suppose: clockwise rotation matrices P (z, −α) 
    R_y = [[cos(rand_y) sin(rand_y) 0];[-sin(rand_y) cos(rand_y) 0];[0 0 1]]; % minus axis Z
%     Ry = roty(ang_y); % Phased Array System Toolbox.% add toolbox of Robot
    y = S(3*i-2:3*i,:);
    S_new = [S_new ; R_y * y ]; 
end

for i=end_pos+1:total_frame
    y = S(3*i-2:3*i,:);
    S_new = [S_new ; R_y * y ];  
end

% plots

x = S(1:3:end, :);
y = S(2:3:end, :);
z = S(3:3:end, :);
xp = [x(:,list(:,1)),x(:,list(:,2))];
yp = [y(:,list(:,1)),y(:,list(:,2))];
zp = [z(:,list(:,1)),z(:,list(:,2))];

% new datas
x_new = S_new(1:3:end, :);
y_new = S_new(2:3:end, :);
z_new = S_new(3:3:end, :);
xp_new = [x_new(:,list(:,1)),x_new(:,list(:,2))];
yp_new = [y_new(:,list(:,1)),y_new(:,list(:,2))];
zp_new = [z_new(:,list(:,1)),z_new(:,list(:,2))];


for i = 1:size(S,1)/3
    subplot(1,2,1)
    title('Original')
    scatter3(x(i,:),y(i,:),z(i,:),'.r');
    for j = 1:length(list)
        hold on;plot3([xp(i,j) xp(i,length(list)+j)],[yp(i,j) yp(i,length(list)+j)],[zp(i,j) zp(i,length(list)+j)],'-',Color='black');
    end
    hold off
    axis equal; drawnow limitrate;
    
    subplot(1,2,2)
    title("With rotations")
    scatter3(x_new(i,:),y_new(i,:),z_new(i,:),'.r');
    for j = 1:length(list)
        hold on;plot3([xp_new(i,j) xp_new(i,length(list)+j)],[yp_new(i,j) yp_new(i,length(list)+j)],[zp_new(i,j) zp_new(i,length(list)+j)],'-',Color='black');
    end
    hold off
    axis equal; drawnow limitrate;
    
end
sgtitle('Data view')
drawnow
view(3);


% for i=1:size(S,1)/3
%     
%     scatter3(S(3*i-2,:),S(3*i-1,:), S(3*i,:),'.r')
%         for j = 1:length(list_p)
%             % add line according to the position of points in the list variable
%             hold on; point_pos = list_p(j,:);
%             p_xyz = [S(3*i-2:3*i,point_pos(1)), S(3*i-2:3*i,point_pos(2))];
%             plot3(p_xyz(1,:),p_xyz(2,:),p_xyz(3,:),'-',Color='black');
%         end
%     hold off
%     axis equal; drawnow;
% %     view(3)
%     
%     subplot(1,2,2)
%     title("With rotations")
%     scatter3(S_new(3*i-2,:),S_new(3*i-1,:), S_new(3*i,:),'filled')
%        for j = 1:length(list_p)
%             % add line according to the position of points in the list variable
%             hold on; point_pos = list_p(j,:);
%             p_xyz = [S_new(3*i-2:3*i,point_pos(1)), S_new(3*i-2:3*i,point_pos(2))];
%             plot3(p_xyz(1,:),p_xyz(2,:),p_xyz(3,:),'-',Color='black');
%         end
%     hold off
%     axis equal; drawnow;
% %     view(3)
%     sgtitle('Data view')
% end
% 


end