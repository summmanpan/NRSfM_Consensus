function S_new = addRotation(ang,S,dataname,list,plot_flag)
% according to the paper:
% starting from the 1/4 position (incluido) of the entire frames
% ending at th2 half position of the frame (excluido)
% the camera rorations were made around the yaxis

total_frame = size(S,1)/3;
start_pos = floor(total_frame/4);% one quarter ej 275.5 = 275
end_pos = start_pos*2; % floor(total_frame/2); <--------------- aqui
ang_per_frame = ang/start_pos; %(end_pos-start_pos-1); 
rad = deg2rad( ang_per_frame );

% rad2 = pi/3; % 60º
% rad = rad2/start_pos;
%inits variables
R_y=0;
rand_y = 0;

% 0 until first quarter frames do not change 
S_new = S(1:(start_pos-1)*3,:); % (start_pos)<-----------------------------aqui

% 1/4 quarter until 2/4 quarter
for i=start_pos:end_pos % start_pos+1:end_pos <--------------aqui
    rand_y = rand_y + rad;
    % R_y = [[cos(rand_y) -sin(rand_y) 0];[sin(rand_y) cos(rand_y) 0];[0 0 1]]; % minus axis Z
    % Suppose: clockwise rotation matrices P (z, −α) 
    R_y = [[cos(rand_y) sin(rand_y) 0];[-sin(rand_y) cos(rand_y) 0];[0 0 1]]; % minus axis Z
    y = S(3*i-2:3*i,:);
    S_new = [S_new ; R_y * y ]; 
end

% remainds frames have last angle
for i=end_pos+1:total_frame %<--------------
    y = S(3*i-2:3*i,:);
    S_new = [S_new ; R_y * y ];  
end

% --plots -> poner en una fuction a parte !?
if ~exist('plot_flag','var')
    return; 
end

% put the access variables separalately to accelerate the coordinate access
x = S(1:3:end, :);
y = S(2:3:end, :);
z = S(3:3:end, :);
x_new = S_new(1:3:end, :);
y_new = S_new(2:3:end, :);
z_new = S_new(3:3:end, :);

if exist('list','var')
    xp_new = [x_new(:,list(:,1)),x_new(:,list(:,2))];
    yp_new = [y_new(:,list(:,1)),y_new(:,list(:,2))];
    zp_new = [z_new(:,list(:,1)),z_new(:,list(:,2))];
    xp = [x(:,list(:,1)),x(:,list(:,2))];
    yp = [y(:,list(:,1)),y(:,list(:,2))];
    zp = [z(:,list(:,1)),z(:,list(:,2))];
end

% plot the original data and with rotation added, per frame in 2 subplots
% if there is the list variable, use it to align the points of each frame
for i = 1:total_frame
    subplot(1,2,1)
    title('Original')
    scatter3(x(i,:),y(i,:),z(i,:),'.r');
    if exist('list','var')
        for j = 1:length(list) % we draw the lines, line by line, according with 2 points
            hold on;plot3([xp(i,j) xp(i,length(list)+j)],[yp(i,j) yp(i,length(list)+j)],[zp(i,j) zp(i,length(list)+j)],'-',Color='black');
        end
        hold off
    end
    axis equal; drawnow limitrate;
    
    subplot(1,2,2)
    title("With rotations")
    scatter3(x_new(i,:),y_new(i,:),z_new(i,:),'filled');
    if exist('list','var')
        for j = 1:length(list) % we draw the lines, line by line, according with 2 points
            hold on;plot3([xp_new(i,j) xp_new(i,length(list)+j)],[yp_new(i,j) yp_new(i,length(list)+j)],[zp_new(i,j) zp_new(i,length(list)+j)],'-',Color='black');
        end
        hold off
    end
    axis equal; drawnow limitrate;
    
end
sgtitle(dataname)
drawnow
view(3);

end