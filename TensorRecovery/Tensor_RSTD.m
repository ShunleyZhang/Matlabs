% T = zeros(5,5,5);
% for i=1:5
% for j=1:5
% for k=1:5
% T(i,j,k)=max(i,j);
% end
% end
% end
% for i= 1:5
%     T(i,6-i,6-i) = 0;
%     T(i,6-i,i) = 0;
% end

  T = double(imread('F:\cadcg410\Desktop\wall.png'));
alpha = [1,1,1e-3];
alpha = alpha / sum(alpha);
beta = alpha;
gamma = alpha*100;
lambda =  alpha*3000;               %lambda up rank down   1000,10000 问题在error的大小和1范数被限制了
eta = 1;
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
stopCriterion = 1;
maxIter = 100;
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
%         disp(temp);
%         disp(tau);
        N_i = max(temp-tau, 0);
        tmps = size(temp);
        for k=1:tmps(1)
            for j=1:tmps(2)
                if (temp(k,j)<-tau)
                    N_i(k,j) = min(temp(k,j)+tau,0);
                end
            end
        end
        
 %       disp(N_i);
        N{i} = N_i;
        %% compute M_i
        L_i = Unfold(L, size(L), i);
        tempuv = (alpha(i)*L_i+gamma(i)*(A{i}-N{i}))/(alpha(i)+gamma(i));
        [U,DIG,V] = svd(tempuv, 'econ');
        
        diagS = diag(DIG);
        tau = lambda(i)/(alpha(i)+gamma(i));
        svp = length(find(diagS > tau));
        sap = length(find(diagS > 0));
  %      fprintf('%d %d\n',svp,sap);
        afterD = max(diagS(1:svp)-tau, 0);
      %  disp(DIG);
       % disp(afterD);
        M_i = U(:,1:svp)*diag(afterD)*V(:,1:svp)';
        %disp(M_i-M{i});
        
%         disp(norm(M_i,  'fro'));
%         disp(norm(M{i},  'fro'));
        %% update
        
        M{i} = M_i;
       % disp(M_i-M{i});
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
% disp(L_star);

subplot(1,3,1);
imshow(uint8(T));
subplot(1,3,2);
imshow(uint8(L_star));
subplot(1,3,3);
imshow(uint8(S_star));
