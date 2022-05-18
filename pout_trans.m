function [X, T] = pout_trans(X, W)
% Eliminate translation component
% Common:
%     X: Trajectory data                              (k x p x f)
%
% Outputs:
%     T: Translation                                  (k x 1 x f)
if nargin < 2
    %X = bsxfun(@minus, X, mean(X, 2));
    T = mean(X, 2);
    X = bsxfun(@minus, X, T);

else
    T = sum(X.*W, 2)./sum(W, 2);
%     X = bsxfun(@minus, X, sum(X.*W, 2)./sum(W, 2)).*W;
    X = bsxfun(@minus, X, T).*W;
    X(isnan(X)) = 0;
end

end


% function [X, T] = pout_trans(X)
% % function [X, T] = pout_trans(X)
% %
% % Eliminate translation component
% %
% % Common:
% %     X: Trajectory data                              (k x p x f)
% %
% % Outputs:
% %     T: Translation                                  (k x 1 x f)
% %
% 
% T = mean(X, 2);
% X = bsxfun(@minus, X, T);
% 
% end