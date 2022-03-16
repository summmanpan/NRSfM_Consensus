function idx = select_idx(D, nsamp, mgroup, lambda)
% function idx = select_idx(D, nsamp, mgroup, lambda, lratio)
%
% Sample trajectory groups
%
% Inputs:
%     D: Input 2D trajectory data                     (k x p x f)
%     nsamp: Number of samples in a trajectory group
%     mgroup: Minimum number of groups that a trajectory belongs to
%     lambda: Parameter for sampling probability
%
% Outputs:
%     idx: Sample indices for trajectory groups       (p x m)
%

    [~, p, f] = size(D);
    
	D = bsxfun(@minus, D, mean(D, 2));

    rD = zeros(size(D));
    for i=1:f
        % para todos los puntos de cada frame
        tD = D(:, :, i); % para todos los puntos 2D de cada frame,
        L = chol(tD*tD'/(p-1), 'lower'); %returns a lower triangular matrix T, such that T*T' = A.
        % chol->Cholesky factorization: para q la matrix inicial se
        % descomponga en two lower dimensionality rectangular matrices
        rD(:, :, i) = L\tD; % matrix cov L dividido invertido entre tD, es decir es realmente
        % tD/L -> entender luego como se div 2Dvector con una matrix!
    end
    rD = reshape(permute(rD, [3 1 2]), [], p); 
    % cambiamos a -> frame x 2 x p , y reshape para que sea N x p , donde
    % p son el numero total de los puntos 
%      frame x 2 x p -> hace que es frame x 2-> es decir, cogemos todos los
%      puntos 2D del primer frame x número de puntos
%   Lo estamos conv en una matrix 2D, donde cada columna contiene la inf de
%   primer frame, donde 1 colum de 1:265 son inf del 1er dimension y luego
%   el de segundo , es decir, [ X' Y' ]

    R = rD'*rD; % estoy representando todos los puntos de la D, factorizando en un simple R

    dR = diag(R);
    R = bsxfun(@plus, bsxfun(@plus, -2*R, dR), dR'); % Suma binaria, pero why??
%     el resultado de esto es que la parte diag son todos 0s.
    R(1:p+1:end) = inf; % cambiamos los elem diag a inf
    lP = R/(-2*nsamp/lambda); % hacemos la normalizacion -> EQ  2!
    lP = lP - max(lP(:)); % log prob->log(P(V_k)) % restamos con el numero max ?whyyy???
    
    idx = select_group(lP, nsamp, mgroup);
    % QUe es idx???
    % Sparse logical, con los coord de los puntos , de 10 puntos, hasta
    % 271?
    %
end

function idx = select_group(lP, nsamp, mgroup)
% Sampling trajectory groups
%
% lP: Log-probability

    p = length(lP);
    
    tID = tic;
    count = 0;
    % S = sparse(i,j,v,m,n,nz)
    idx = sparse(zeros(1, 0), zeros(1, 0), false(1, 0), p, mgroup*p, nsamp*mgroup*p); %???
    % idx = sparse ([], [], empty.logical, 91, 50*91=4550, 10*50*91=45500
    % ) = 91*4550 sparse logical

    % []-> empty double row vector
    % empty logical array
    %converts a full matrix into sparse form by squeezing out any zero elements. 
    % If a matrix contains many zeros, converting the matrix to sparse storage saves memory.

    npidx = zeros(p, 1); % numb of points 
    while any(npidx < mgroup) % si hay algun numero menor que mgrupo
        % si todos son mayores que mgrup -> break

        % supongo es para asegurar de que todos los puntos are included in
        % at least m_g????
        count = count + 1;
        
        tlP = zeros(1, p);
        j = random_select(log(double(npidx == min(npidx)))); % log(double(npidx == min(npidx)))?
        idx(j, count) = true;
        npidx(j) = npidx(j)+1;
        tgind = idx(j, 1:count-1); % coge valores de idx, de la fila j todas las columnas hasta count
        tlP = tlP + lP(j, :); % esta acumulando las probabilidades, de gasussian
        % first traj is based on a uniform distribution
        % ahora, second traj based on gaussian!
        for i=2:nsamp 
            tmp = tlP > -inf; %comprobamos q sea mayores que -inf
            % why comprobamos que una prob sea mayor que -inf???
            % HE REPASADO HASTA AQUI!!!!!!!
            tmp(tmp) = check_comb( sum(idx(tmp, tgind), 2) , p-i, nsamp-i);
            tlP(tmp) = -inf; % tranf a inf los elementos con index tmp,
            %como no hay posiscion 0,0, para la primera no hay -inf ?
            
            [j, tlP] = random_select(tlP);
            idx(j, count) = true;
            npidx(j) = npidx(j)+1;
            tgind = tgind & idx(j, 1:count-1);
            tlP = tlP + lP(j, :);
        end
        
        if toc(tID) > 1
            disp(['select_idx ' num2str(count) ' / ' num2str(min(npidx)) ' / ' num2str(mean(npidx)) ' / ' num2str(mgroup)]);
            tID = tic;
        end
    end
    
    idx = idx(:, 1:count); %??
end


function [i, lP] = random_select(lP)
% Sampling a trajectory
%
% i: Index of the sampled trajectory

    lP = lP - max(lP);
    cP = cumsum(exp(lP)); % weight sampling???
    lP = lP - log(cP(end)); % coge el ultimo num, hace el log, y lo resta.
    i = sum(cP < rand*cP(end)) + 1;
end


function r = check_comb(x, n, k)
% Rule out unnecessary combinations
%
% r = (x == nchoosek(n, k))
% r: 1 if there are no remaining combinations
%
% Check whether there are no remaining combinations to speed up the process

    if log(max(x)) > -betaln(n-k+1, k+1) - log(n+1)-1 
    % If the log of maximum is larger than log(nchoosek/exp) 
    % log(nchoosek/exp) ?????
        r = (x == nchoosek(n, k));
    %nchoosek(N,K) where N and K are non-negative integers returns N!/K!(N-K)!.
%     is the number of combinations of N things taken K at a time.
    else
        r = false(size(x));
    end
end




