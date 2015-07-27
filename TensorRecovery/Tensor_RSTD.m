T = zeros(5,5,5);
for i=1:5
for j=1:5
for k=1:5
T(i,j,k)=max(i,j);
end
end
end
T(1,2,3) = 0;
% T = double(imread('E:\Users\Desktop\Tensor\testImg.png'));
alpha = [1,1,1];
beta = alpha;
gamma = alpha/10;
lambda =  alpha;
eta = 0.1;
%T -    n mode tensor of observations/data (required input)
%alpha
%beta
%gamma
%lambda
%eta

%% initialization
A = cell(ndims(T), 1);
M = cell(ndims(T), 1);
N = cell(ndims(T), 1);
L = T;
sz = size(T);
S = T;  
S(:,:,:) = 0;%for 3d tensor
for i = 1:ndims(T)
        dim = size(T);
        M{i} = Unfold(T, dim, i);
        A{i} = Unfold(T, dim, i);
        sz = size(M);
        N{i} = Unfold(S, dim, i);
end
k = 1;
t = 1;

iter = 0;
converged = false;
stopCriterion = 3;
maxIter = 3;
tol = 1e-14;
Tnorm = norm(Unfold(T,size(T),1), 'fro');
%% main loop

while ~converged       
    iter = iter + 1;
    
    for i = 1:ndims(T)
        %% compute N_i
        S_i = Unfold(S, size(S), i);
        temp = (beta(i)*S_i+gamma(i)*(A{i}-M{i}))/(beta(i)+gamma(i));
        tau = eta/(beta(i)+gamma(i));
        N_i = max(temp-tau, 0);
        %% compute M_i
        L_i = Unfold(L, size(L), i);
        tempuv = (alpha(i)*L_i+gamma(i)*(A{i}-N{i}))/(alpha(i)+gamma(i));
        [U,DIG,V] = svd(tempuv, 'econ');
      % disp(DIG);
        diagS = diag(DIG);
        tau = lambda(i)/(alpha(i)+gamma(i));
        svp = length(find(diagS > tau));
        sap = length(find(diagS > 0));
        fprintf('%d %d\n',svp,sap);
        afterD = max(diagS(1:svp)-tau, 0);
        M_i = U(:,1:svp)*diag(afterD)*V(:,1:svp)';
        %disp(M_i);
        
        disp(norm(M_i,  'fro'));
        disp(norm(M{i},  'fro'));
        %% update
        N{i} = N_i;
        M{i} = M_i;
    end
    
    %% compute L and S
    S_star = beta(1) * Fold(N{1}, size(T), 1);
    S_sum = beta(1);
    L_star = alpha(1) * Fold(M{1},size(T), 1);
    L_sum = alpha(1);
    for i = 2:ndims(T)
        S_star = S_star + beta(i) * Fold(N{i},size(T), i);
        S_sum = S_sum + beta(i);
        L_star = L_star + alpha(i) * Fold(M{i},size(T), i);
        L_sum = L_sum + alpha(i);
    end
    
    %% stop Criterion 
    Error = T - S_star - L_star;
    stopCriterion = norm(Unfold(Error, size(T), 1), 'fro') / Tnorm;
    if stopCriterion < tol
        converged = true;
    end    
    
     disp(['Iteration' num2str(iter) ' stopCriterion ' num2str(stopCriterion)]);
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

