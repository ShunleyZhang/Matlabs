

%T -    3 mode tensor of observations/data (required input)
%alpha
%beta
%gamma
%lambda
%eta
%rho
%L
%S

%% initialization
rho = 1;
k = 1;
L = cell(ndims(T), 1);
M = cell(ndims(T), 1);
N = cell(ndims(T), 1);
for i = 1:ndims(T)
        L{i} = Unfold(T, dim, i);
        M{i} = L{i};
        sz = size(M{i});
        N{i} = zeros(sz(1), sz(2));
end
%% main loop
while ~converged       
    for i = 1:ndims(T)
        temp = (alpha(i)*(L{i}+Y{i}/alpha(i))
    end
    
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = sv + round(n * 0.1);
    while primal_converged == false
        
        temp_T = D - A_hat + (1/mu)*Y;
        temp_E = max( temp_T - lambda/mu,0) + min( temp_T + lambda/mu,0); 
        
        if choosvd(n, sv) == 1
            [U S V] = lansvd(D - temp_E + (1/mu)*Y, sv, 'L');
        else
            [U S V] = svd(D - temp_E + (1/mu)*Y, 'econ');
        end
        diagS = diag(S);
        svp = length(find(diagS > 1/mu));
        if svp < sv
            sv = min(svp + 1, n);
        else
            sv = min(svp + round(0.05*n), n);
        end
        temp_A = U(:,1:svp)*diag(diagS(1:svp)-1/mu)*V(:,1:svp)';    
        
        if norm(A_hat - temp_A, 'fro') < tolProj && norm(E_hat - temp_E, 'fro') < tolProj
            primal_converged = true;
        end
        A_hat = temp_A;
        E_hat = temp_E;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
   
    %**********************
    S_star =1;
    L_star =1;
    Error = T - S_star - L_star;
    
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / dnorm;
    if stopCriterion < tol
        converged = true;
    end    
    
    disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' r(A) ' num2str(svp)...
        ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
        ' stopCriterion ' num2str(stopCriterion)]);
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end