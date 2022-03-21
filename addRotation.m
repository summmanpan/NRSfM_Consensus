function GT = addRotation(ang,GT,dataname)
% according to the paper
% starting from the 1/4 position of the entire frames
% ending at th2 half position of the frame
% the camera rorations were made around the yaxis

% syms R_y(ang)
% R_y(ang) = [[cos(ang) 0 sin(ang)];[0 1 0];[-sin(ang) 0 cos(ang)]];

total_frame = size(GT,3);
start_pos = round(total_frame/4);
end_pos = round(total_frame/2);

ang_per_frame = ang/total_frame; % todo el frame
% ang_per_frame = ang/(end_pos-start_pos); % solo el cuarto parte

ang_y = 0;
rand_y = 0;
rad = deg2rad( ang_per_frame );
GT_R_manual = GT;

for i=1:total_frame %i=start_pos:end_pos % per frame
    rand_y = rand_y + rad;
%     ang_y = ang_y + ang_per_frame;
    % Ry = R_y(ang_per_frame);
    R_y = [[cos(rand_y) 0 sin(rand_y)];[0 1 0];[-sin(rand_y) 0 cos(rand_y)]];

%     Ry = roty(ang_y); % Phased Array System Toolbox.% add toolbox of Robot
    

%     for each points
%     for j=1:size(frame,2)
%         y = frame(1,j,i); 
%         quat = quaternion([0,ang_per_frame,0],'eulerd','XYZ','frame');
%         rereferencedPoint_y = rotateframe(quat,[0 y 0]);
% %         guardar en , si quiero ya lo miro si vale la pena o no hacerlo con cuaternion!
%         
%     end


    y = GT(1:3,:,i);
    GT(1:3,:,i) = R_y * y; 
    GT_R_manual(1:3,:,i) = inv(R_y)*y ; 

    % rereferencedPoint = rotateframe(quat,[x,y,z])
    
end



for i=1:total_frame
    scatter(GT(1, :, i), GT(2, :, i) , 'ro');
    hold on; 
    scatter(GT_R_manual(1, :, i), GT_R_manual(2, :, i) , 'bo','filled');
    hold off;
    axis equal; title(dataname); drawnow;
end             

% serch ang btw 2 vectors
% CosTheta = max(min(dot(vec1,vec2)/(norm(vec1)*norm(vec2)),1),-1);
% ThetaInDegrees = real(acosd(CosTheta))

end