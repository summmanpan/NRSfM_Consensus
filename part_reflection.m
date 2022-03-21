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
    % cada Zi o Z_grupo o Zk guarda el 3º componente en col de todos los
    % frames para cada grupo i, X{i}
end
Z = reshape(Z, [], n);

% hacemos la normalizacion con su modulo
Z = bsxfun(@rdivide, Z, sqrt(sum(Z.^2))); 
% Z / ||Z||, y por q hacemos esto?
%sum(X) is the sum of the elements of the vector X. If X is a matrix,
%     S is a row vector with the sum over each column
% sum(,2) es la suma por col, ie. queda una col con todas las sumas de cada
% vector fila

nmax = max(sum(idx, 2));
% esto es la suma de para cada uno de los puntos p, el punto que ha
% parecido mas veces en el random sampling¿¿??
% y vemos q al menos el nº de apariciones son mayor que 50.
% que significa el número de ....???

if min(size(Z)) < 1500 % compara el n < 1500
    [U, ~, ~] = svd(Z, 'econ');
    U = U(:, 1);
    %produces the "economy size" decomposition.
    % if Z is this case is m>=n, then ç
    %only the first n columns of U are computed and 
    % S is n-by-n.

else
    %Find a few singular values and vectors.
    %svds(A,K):computes the K largest singular values of A.
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

