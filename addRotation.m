function S_new = addRotation(ang,S,dataname)
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

for i=1:size(S,1)/3
%     plot3(S(3*i-2,:),S(3*i-1,:),S(3*i,:),'.r');
    subplot(1,2,1)
    title("Original")
%     scatter3(S(3*i-2,:),S(3*i-1,:), S(3*i,:),'.r')
    plot3(S(3*i-2,:),S(3*i-1,:),S(3*i,:),'.r');
    view(3)
    axis equal; drawnow;

    subplot(1,2,2)
    title("With rotations")
    scatter3(S_new(3*i-2,:),S_new(3*i-1,:), S_new(3*i,:),'filled')
    view(3)
    axis equal; drawnow;
end



end