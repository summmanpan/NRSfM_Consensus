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
    
	D = bsxfun(@minus, D, mean(D, 2)); % eliminate translation
    % D = D - mean(D,2); (they have slightly -4.1599e-12 differences)

    rD = zeros(size(D));
    for i=1:f
        
        tD = D(:, :, i); % para todos los puntos 2D de cada frame,
        L = chol(tD*tD'/(p-1), 'lower'); %returns a lower triangular matrix T, such that T*T' = A.
        % chol->Cholesky factorization: para q la matrix inicial se
        % descomponga en two lower dimensionality rectangular matrices
        rD(:, :, i) = L\tD; %  mldivide(A,B)
        % A\B is the matrix division of A into B,
% %     X = A\B is the solution to the equation A*X = B.
    end
    rD = reshape(permute(rD, [3 1 2]), [], p); 
    % cambiamos -> frame x 2 x p , y reshape para que sea N x p , donde
    % p son el numero total de los puntos 
    %el permute hace que guardemos por puntos, como varian los puntos a lo
    %largo del todo el frame. 1º columna inf de Xs, y luego Ys, [X', Y']
    % el reshape hace que guardemos por col,toda la inf vectorizado del 1
    % pto [X' Y']

    R = rD'*rD; % estoy representando todos los puntos de la D, factorizando en un simple R!

    dR = diag(R);
    R = bsxfun(@plus, bsxfun(@plus, -2*R, dR), dR');% suma por filas 
    %bsxfun(@plus, -2*R, dR)-> 1 fila de dR suma a 1 fila de todos -2R ->
    %la parte diag son todos igual que antes, cambiando el signo
%     el resultado final de esto es que la parte diag son todos 0s.
    R(1:p+1:end) = inf; % cambiamos los elem diag a inf
    lP = R/(-2*nsamp/lambda); % hacemos la normalizacion -> EQ  2!
    lP = lP - max(lP(:)); % log(P(V_k)) % restamos con el numero max 
    
    idx = select_group(lP, nsamp, mgroup);
    %  idx es el sparse matrix donde guardo todos mis clusters!
    % 55 son el número total de puntos para mi imagen...
    % Sparse logical, con los coord de los puntos , de 10 puntos, hasta
    % 55x2750 hasta ser 55x362 columnas
    % contiene los 362 grupos, para cada 10 trajectorias.
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
    % idx = sparse ([], [], empty.logical, 55, 50*55, 10*50*55) = 55*2750 sparse logical

    % []-> empty double row vector empty logical array
    % converts a full matrix into sparse form by squeezing out any zero elements. 
    % If a matrix contains many zeros, converting the matrix to sparse storage saves memory.

    npidx = zeros(p, 1); % numb of points 
    while any(npidx < mgroup) % si hay algun numero menor que mgrupo
        % si todos son mayores que mgrup -> break

        % supongo es para asegurar de que todos los puntos are included in
        % at least m_g???? Es decir, el punto sampleado tiene que estar
        % repetido al menos en 55 grupos?

        count = count + 1;
        
        tlP = zeros(1, p);
        j = random_select(log(double(npidx == min(npidx)))); %idx que no congentan los idx encontrados anterirmente
        idx(j, count) = true;
        npidx(j) = npidx(j)+1; % el vector con indx pos flags, necesitamos que cada pto tenga al menos mgroup muestras.
        tgind = idx(j, 1:count-1); %de idx, coge de la fila j todas las columnas hasta count-1. el punto p, los grupos que perteneció!
        tlP = tlP + lP(j, :); % esta acumulando las probabilidades, lo guarda actualizando para cada p,y justo la fila j, tendrá -inf
        % first traj is based on a uniform distribution
        
        % ahora, second traj based on gaussian!
        % 
        for i=2:nsamp % para los 10 samples de este grupo?
            tmp = tlP > -inf; %comprobamos q sea mayor que -inf
            % asi comprobamos que la prb no sera muy peq?
            
            tmp(tmp) = check_comb( sum(idx(tmp, tgind), 2) , p-i, nsamp-i);%idx(tmp, tgind) busca las filas que sean empty tgind
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
    
    idx = idx(:, 1:count); 
end


function [i, lP] = random_select(lP)
% Sampling a trajectory
%
% i: Index of the sampled trajectory

    lP = lP - max(lP);
    cP = cumsum(exp(lP)); % vector of (1 to 55)
    % ->weight sampling???
    % for the fisrt one. since lP is 0, e^0 es 1. therefore cP cumsum of
    % 1:55
    lP = lP - log(cP(end)); % coge el ultimo num, hace el log, y lo resta.
    i = sum(cP < rand*cP(end)) + 1; % 
    % rand -> uniform distribution [0 1]*55, --> vales btw [0 55]
    %sum(cP < rand*cP(end)) compara con el vector, y lo sumamos y +1 ???
    % why we do that??
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




