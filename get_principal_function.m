function [X, err] = get_principal_function(GT,dataname,flag_regu, regu_type, regu_order,noise)
% GT: data with gt of 3d structure
% dataname: sting of filename 
% flag_regu: FLAG OF APPLY REGU OR NOT --- 1,2 INTEGER
% regu_type: REGU TYPE ---SOFT OR HARD --- string
% regu_order: ORDER OF REGULARIZATION--- 1,2,4---INTEGER

% noise, miss data generation.
% generate the diary and calculate the error

%% Input data generation

% Experimental setting
% noise = 0; %10^-3;            % noise level. paper use 10^-3

[k, p, nSample] = size(GT);
D = zeros(k, p, nSample);
temp = GT(1:2, :, :);

% add normal random number
weight_noise = noise*max(abs(reshape(bsxfun(@minus, temp, mean(temp, 2)), [], 1)));
D(1:2, :, :) = temp + weight_noise*randn(2, p, nSample);

% Consensus of Non-Rigid Reconstructions
X = NRSfM_Consensus(D, flag_regu, regu_type, regu_order );

% save("back_sparse_reconst_with_Regu_soft_l2",'X')
err=0;

if dataname == "back_sparse"
    save("./back_plots/back_sparse_reconst_with_C_Hard_2",'X')
    return
end
%% Evaluation Error

GT = bsxfun(@minus, GT, mean(GT, 2));
vind = sum((GT(3, :, :)-X(3, :, :)).^2) > sum((GT(3, :, :)+X(3, :, :)).^2);
% > q hay mas diff GT que X, es decir si GT es 10, y X 7, pues la diff si
% es mas mas grande que la suma de 10 mas 7 , quiero decir que
% si es 10-10=0 > 20-> no errror
% si error es mayor que 0, o no error, lo invertimos.
% si es 3 > ?? 1-10=-9 > 10 -> falso-> es neg -> invertimos
% aseguramos de que la diff de GT - X, es positivo 

X(3, :, vind) = -X(3, :, vind);

perf = sqrt(reshape(sum(sum((GT-X).^2)), 1, [])./reshape(sum(sum(GT.^2)), 1, []));

err = num2str(mean(perf)); 
disp(['---' dataname '---------MEAN ERROR---------'])
disp(['mean error : ' err]);



end