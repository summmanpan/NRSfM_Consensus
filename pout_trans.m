function [X, T] = pout_trans(X)
% function [X, T] = pout_trans(X)
%
% Eliminate translation component
%
% Common:
%     X: Trajectory data                              (k x p x f)
%
% Outputs:
%     T: Translation                                  (k x 1 x f)
%

T = mean(X, 2);
X = bsxfun(@minus, X, T);

end