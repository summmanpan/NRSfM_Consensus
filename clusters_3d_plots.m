function clusters_3d_plots(Xi)
% Plot Xi clusters results before the consensus

% Inputs:
%     Xi: all the clusters  (1 x f) cell

hFig = figure();
axh = axes('Parent', hFig);
hold(axh, 'all');
grid(axh, 'on');
frames_number = size(Xi{1}(:,:,:),3);
for i=1: frames_number
    number_end =15;
    for s = 1:2
        if s == 1
            subplot(1,3,1)
        else
            subplot(1,3,2)
        end
        for j = 1:number_end %
            r_part = Xi{j}(:,:,:); 
            scatter3(r_part(1, :, i), r_part(3, :, i), r_part(2, :, i),'filled');
            hold on
            if j==number_end % llega al final
                hold off
            end
        end
       
        if s == 1
           %subplot(1,3,3)
%            hold on
%            scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
           view(3);
%            axis equal; drawnow;
           %view([AZ,EL])
           hold off
        else
            view(2); 
        end
        axis equal; drawnow;
    end
%     subplot(1,3,3)
%     scatter3(GT(1, :, i), GT(3, :, i), GT(2, :, i), 'ro');
%     view(3);
%     axis equal; drawnow;  
end

xlabel('X')
ylabel('Y')
zlabel('Z')

end