function plot_NRSfM(D, W, gX, rX, vidObj)
% function plot_NRSfM(D, W, gX, rX, vidObj)
%
% Plot reconstructed results
%
% Inputs:
%     D: Observations                                 (3 x p x n)
%     W: Weights - observed (true) or missing (false) (3 x p x n)
%     gX: Ground truth shapes                         (3 x p x n)
%     rX: Reconstructed shapes                        (3 x p x n)
%     vidObj: VideoWriter object (write a video file if specified)
%
% Ref: Minsik Lee, Jungchan Cho, Chong-Ho Choi, and Songhwai Oh,
% "Procrustean Normal Distribution for Non-Rigid Structure from Motion,"
% CVPR 2013, Portland, Oregon, June 23-28, 2013.
%
% Author: Minsik Lee (mlee.paper@gmail.com)
% Last update: 2013-09-07
% License: GPLv3

%
% Copyright (C) 2013 Minsik Lee, Jungchan Cho
% This file is part of NRSfM_PND.
%
% NRSfM_PND is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% NRSfM_PND is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with NRSfM_PND.  If not, see <http://www.gnu.org/licenses/>.


if nargin < 5 || isempty(vidObj)
    f_rate = 30;
    v_flag = false;
else
    open(vidObj) %
    f_rate = vidObj.FrameRate;
    v_flag = true;
end

nSample = size(D, 3);

D = bsxfun(@minus, D, sum(D.*W, 2)./sum(W, 2)).*W;
D(isnan(D)) = 0;
gX = bsxfun(@minus, gX, mean(gX, 2));
rX = bsxfun(@minus, rX, mean(rX, 2));
ind = sum((gX(3, :, :)-rX(3, :, :)).^2) > sum((gX(3, :, :)+rX(3, :, :)).^2);
rX(3, :, ind) = -rX(3, :, ind);

axD = [min(D(1, W(1, :))) max(D(1, W(1, :))) min(D(2, W(2, :))) max(D(2, W(2, :)))];
axX = [min([gX(1, :) rX(1, :)]) max([gX(1, :) rX(1, :)]) min([gX(3, :) rX(3, :)]) max([gX(3, :) rX(3, :)]) min([gX(2, :) rX(2, :)]) max([gX(2, :) rX(2, :)])];

h = clf('reset');
set(h, 'Color', 'w');

for k=1:nSample
    TT = tic;
    set(h, 'Name', [num2str(k) ' / ' num2str(nSample)]);
    
    subplot('Position',[0 0 1/3 0.9]);
    ind = all(W(1:2, :, k));
    scatter(D(1, ind, k), D(2, ind, k), 'k.');
    axis('equal', axD, 'off'); grid off;
    title('Input 2D tracks');

    subplot('Position',[1/3 0 1/3 0.9]);
    scatter3(gX(1, :, k), gX(3, :, k), gX(2, :, k), 'ro'); hold on;
    scatter3(rX(1, :, k), rX(3, :, k), rX(2, :, k), 'b+'); hold off;
    axis('equal', axX, 'off'); grid off;
    view(45, 30); title('3D view 1');
    
    subplot('Position',[2/3 0 1/3 0.9]);
    scatter3(gX(1, :, k), gX(3, :, k), gX(2, :, k), 'ro'); hold on;
    scatter3(rX(1, :, k), rX(3, :, k), rX(2, :, k), 'b+'); hold off;
    axis('equal', axX, 'off'); grid off;
    view(-45, 30); title('3D view 2');
    
    legend('Ground truth', 'Reconstructed', 'Location', 'SouthEast');
    drawnow;

    if v_flag
        writeVideo(vidObj, getframe(h));
    end
    
    pause(1/f_rate-toc(TT));
end

%
if v_flag
    close(vidObj);
end

end


