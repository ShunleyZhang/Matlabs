%%Test single channel demage 
I = imread('fruit.jpg');

T = I;

y = wgn(size(I,1),size(I,2),0.2);
y = max(y,-1);
y = min(y,1);

for i = 1:size(I,1)
    for j = 1:size(I,2)
        if abs(y(i,j))<1
            T(i,j,1) = min(T(i,j,1) * (1+y(i,j)),255);
        end
    end
end

R = (I>=0);
Omega0 = (abs(y)>=1);
R(:,:,1) = Omega0;
dT = double(T);
dI = double(I);
%[A_hat,E_hat] = exact_alm_rpca(doubleC, 0.05);
alpha = [1, 1, 1e-3];
alpha = alpha / sum(alpha);
maxIter = 500;
epsilon = 1e-5;
mu = 5 * alpha ./ sqrt(size(dT));
C =  0.6;
L0 = 1e-5;
[X_F, errList_F] = FaLRTC(...
     dT,...                    % a tensor whose elements in Omega are used for estimating missing value
     R,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
     );
 
subplot(221);
imshow(I);
subplot(222);
imshow(T);
subplot(223);
imshow(T(:,:,1));
subplot(224);
imshow(uint8(X_F));
