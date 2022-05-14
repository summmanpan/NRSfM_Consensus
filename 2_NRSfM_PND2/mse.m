function result = mse(X)
% Calculate mean squared error
% (For those who do not have Neural Network Toolbox)
% 
% We would like to thank Prof. Chan Su Lee for the valuable comment on the compatibility of the code.
%
% Implemented by Minsik Lee (mlee.paper@gmail.com)
% Last update: 2013-08-15

    result = mean(X(:).^2);
end
