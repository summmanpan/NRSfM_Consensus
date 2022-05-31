function plot_2D(X, typedata)
% Plot 2D information
% X: is the 3D information
disp("The options of the plots are: ")
disp('1: scatter, 2: 1by1, 3: allin1')
n = input('Enter a number of plot you want: ');
X = X(1:2, :, :);

% cellImgStore = {};
for i=1:size(X,3)
    switch n
        case 1
            
            scatter(X(1, :, i), X(2, :, i), 25, ...
            'MarkerEdgeColor',[0 .5 .5],...
            'MarkerFaceColor',[0 .7 .7],...
            'LineWidth',1.5);
            axis equal; 
            title(typedata);
            drawnow;
        case 2

            clf;
            scatter(X(1, :, i), X(2, :, i), 25, ...
                    'MarkerEdgeColor',[0 .5 .5],...
                    'MarkerFaceColor',[0 .7 .7],...
                    'LineWidth',1.5);
            hold on;
            title(typedata);
            getframe;
            pause;
        case 3
%             clf;
%            %# show them in subplots
%             figure(1)
%             for ii=1:6
%                 subplot(2,3,ii);
%                 img = scatter(X(1, :, ii), X(2, :, ii), 25, ...
%                 'MarkerEdgeColor',[0 .5 .5],...
%                 'MarkerFaceColor',[0 .7 .7],...
%                 'LineWidth',1.5);axis equal; 
% %                 cellImgStore{end+1}= figure1;
%                 set(img, 'ButtonDownFcn',{@callback,i})
%             end
        otherwise
            
    end     
end



end
