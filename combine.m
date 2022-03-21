function [X, T] = combine(D, Xi, idx)
% function [X, T] = combine(D, Xi, idx)
% Obtain strong reconstruction
% Inputs:
%     D: Input 2D trajectory data                     (k x p x f)
%     Xi: 3D trajectory groups (corrected sign)       (1 x m cell)
%     idx: Sample indices for trajectory groups       (p x m)
% Outputs:
%     X: Strong 3D reconstruction                     (k x p x f)
%     T: Corrected signs of trajectory groups         (m x f)

    f = size(Xi{1}, 3);
    p = size(idx, 1);
    
    tID = tic;
    Z = zeros(p, f);
    T = zeros(numel(Xi), f);
    [Az, At] = build_A(idx); %Az??
    A = [Az At];
    Ap = A(1:end-1, :);
    AAp = Ap'*Ap;
    lp = eigs(A'*A, 1); % el maximo eigenvalue
    for l=1:f
        Yl = build_Y(Xi, l);
        x = solve(Yl, A, AAp, lp, p); % ADMM problem
        Z(:, l) = x(1:p);
        T(:, l) = x(p+1:end);
        
        if toc(tID) > 1 % if Elapsed time  pasa mas de 1 min , pero porque asi??
            disp(['combine ' num2str(l) ' / ' num2str(f)]);
            % ha combinado 1 de 260 frames
            tID = tic;
        end
    end
    
    X = D;
    X(size(D, 1), :, :) = Z;

end

function x = solve(y, A, AAp, lp, p)
% Solve ADMM problem
%
% minimize || y - A x ||_1
% AAp: Part of A'*A (precomputed for efficiency)
% lp: Lipschitz constant
% p: Number of landmarks
% x: Reconstructed frame (Z and T)

    mu = 1e-4;
    rho = 1.01;
    
    x = zeros(length(AAp), 1);
    L = zeros(numel(y), 1);
    yAx = y;
    
    count = 0;
    cost = inf;
    while cost > 1e-10
        
        tmp = yAx + L;
        E = sign(tmp).*max(abs(tmp)-1/mu, 0);
        
        AAx = AAp*x;
        AAx(1:p) = AAx(1:p) + sum(x(1:p));
        x = x - (AAx - A'*(y - E + L))/lp;
        yAx = y-A*x;
        
        Q = yAx - E;
        L = L + Q;
        cost = norm(Q(:))^2;

        L = L/rho;
        mu = mu*rho;
        
        count = count + 1;
%         if mod(count, 1e2) == 0
%             disp(['solve ' num2str(count) ' : ' num2str(sum(abs(E(:)))) ' / ' num2str(cost) ' / ' num2str(mu)]);
%         end
    end
%     disp(['solve ' num2str(count) ' : ' num2str(sum(abs(E(:)))) ' / ' num2str(cost) ' / ' num2str(mu)]);
    
end

function Y = build_Y(Xi, l)
% Concatenate the reconstructions of different groups
%
% l: Frame index
% Y: Concatenated data

    d = size(Xi{1}, 1); % d es la dimension, en este caso la z, el 3er axis
    Y = cell2mat( cellfun(@(x) { x(d, :, l)'}, reshape(Xi, [], 1)) );
%     Xi = 1x361, Y = 361*10 = 3610x1 size, ya q hay 10 pts
    Y = [Y; 0];
end

function [Az, At] = build_A(idx)
% Create A matrix
%
% Az: Part for z-coordinates
% At: Part for translations

    Az = cell2mat(cellfun(@(x) {selection_A(x)}, num2cell(idx, 1)'));
    %A = cellfun(FUN, C) applies the function specified by FUN to the
%     contents of each cell of cell array C
    tmp = arrayfun(@(n) {sparse(1:n, 1, 1)}, full(sum(idx)));
    At = blkdiag(tmp{:});
    
    Az = [Az; sparse(1, 1:size(Az, 2), 1)];
    At = [At; sparse(1, size(At, 2))];
end

function E = selection_A(idx)
% Create a mapping matrix
    s = sum(idx); %???
    E = sparse(1:s, find(idx), 1, s, numel(idx));
    
end