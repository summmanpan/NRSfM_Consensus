function r = part_reflection(X, idx)
% function r = part_reflection(X, idx)
%
% Resolve reflection ambiguities between groups
%
% Inputs:
%     X: 3D trajectory groups                         (1 x m cell)
%     idx: sample indices for trajectory groups       (p x m)
%
% Outputs:
%     r: corrected signs of trajectory groups         (1 x m)
%

f = size(X{1}, 3);
[p, n] = size(idx);
nidx = sum(idx);
Wp = zeros(p);
for i=1:n
    Wp(idx(:, i), idx(:, i)) = Wp(idx(:, i), idx(:, i)) + eye(nidx(i)) - ones(nidx(i))/nidx(i);
end

Z = zeros(p, f, n);
for i=1:n
    Z(idx(:, i), :, i) = X{i}(3, :, :);
end
Z = reshape(Z, [], n);
Z = bsxfun(@rdivide, Z, sqrt(sum(Z.^2)));

nmax = max(sum(idx, 2));
if min(size(Z)) < 1500
    [U, ~, ~] = svd(Z, 'econ');
    U = U(:, 1);
else
    [U, ~, ~] = svds(Z, 1);
end

tID = tic;
pU = zeros(size(U));
while mse(pU(:)-U(:)) > eps
    pU = U;
    
    U = nmax*U + Z*(U'*Z)' - reshape(Wp*reshape(U, p, f), [], 1);
    U = U/norm(U);
    
    if toc(tID) > 1
        disp(['part_reflection ' num2str(mse(pU-U))]);
        tID = tic;
    end
end

r = (U'*Z >= 0)*2-1;

end

