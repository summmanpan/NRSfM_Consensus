function GT = addRotation(ang,GT)

% syms R_y(ang)
% R_y(ang) = [[cos(ang) 0 sin(ang)];[0 1 0];[-sin(ang) 0 cos(ang)]];



% according to the paper
% starting from the 1/4 position of the entire frames
% ending at th2 half position of the frame
% the camera rorations were made around the yaxis

total_frame = size(GT,3);
start_pos = round(total_frame/4);
end_pos = round(total_frame/2);

ang_per_frame = ang/(end_pos-start_pos);
% quat = quaternion([0,ang_per_frame,0],'eulerd','XYZ','frame');
% R = roty(ang_per_frame); % Phased Array System Toolbox.
% Ry = R_y(ang_per_frame);

ang_y = 0;
% Y_rotated = zeros(size(GT));
% GT(3,:,:) = 0; % eliminate z-axis %?????? elimino?
rad = deg2rad( ang_per_frame );
for i=1:total_frame %i=start_pos:end_pos % per frame
    ang_y = ang_y + rad;
    R_y = [[cos(ang_y) 0 sin(ang_y)];[0 1 0];[-sin(ang_y) 0 cos(ang_y)]];

    % for each points
%     for j=1:size(frame,2)
%         %x = frame(1,j,i);
%         y = frame(1,j,i);
% %         rereferencedPoint = rotateframe(quat,[x,y,z]);
%         % rotation matrix R and vector y
%         frame_rotatedY() = Ry * y;
%     end
    y = GT(1:3,:,i);
    GT(1:3,:,i) = inv(R_y) * y; % luego comprobar esto si realmente esta bien !!!!
    
end


% R = rotx(ang);
% Y_rotated = R*Y;
% add toolbox of Robot...
               
% rereferencedPoint = rotateframe(quat,[x,y,z])
% 
% plot(rereferencedPoint(1,1),rereferencedPoint(1,2),'bo')
% plot(rereferencedPoint(2,1),rereferencedPoint(2,2),'go')

end