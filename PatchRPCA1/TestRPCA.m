%%Test single channel demage 
I = imread('testImg.jpg');
I1 = im2double(I);
dI = I1(:,:,1);
[A_hat,E_hat] = exact_alm_rpca(dI, 0.05);
imshow(A_hat);
% 
% T = I;
% 
% y = wgn(size(I,1),size(I,2),1);
% y = max(y,-1);
% y = min(y,1);
% 
% tmp = 0;
% for i = 1:size(I,1)
%     for j = 1:size(I,2)
%         if abs(y(i,j))<1
%             T(i,j,1) = min(T(i,j,1) * (1+y(i,j)*5),255);
%             tmp = tmp+1;
%         end
%     end
% end
% disp(tmp/size(I,1)/size(I,2));
% R = (I>=0);
% Omega0 = (abs(y)>=1);
% R(:,:,1) = Omega0;
% dT = double(T);
% dI = double(I);

%[A_hat,E_hat] = exact_alm_rpca(doubleC, 0.05);

%%test tensor completion
% alpha = [1, 1, 1e-3];
% alpha = alpha / sum(alpha);
% maxIter = 500;
% epsilon = 1e-5;
% mu = 5 * alpha ./ sqrt(size(dT));
% C =  0.6;
% L0 = 1e-5;
% [X_F, errList_F] = FaLRTC(...
%      dT,...                    % a tensor whose elements in Omega are used for estimating missing value
%      R,...          % the index set indicating the obeserved elements
%      alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
%      mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
%      L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
%      C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
%      maxIter,...        % the maximum iterations
%      epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
%      );
%  
% subplot(241);
% imshow(I);
% subplot(242);
% imshow(T);
% subplot(243);
% imshow(uint8(X_F));
% subplot(244);
% imshow(Omega0);
% 
% subplot(245);
% imshow(I(:,:,1));
% subplot(246);
% imshow(T(:,:,1));
% subplot(247);
% imshow(uint8(X_F(:,:,1)));
% subplot(248);
% imshow(abs(I(:,:,1)-uint8(X_F(:,:,1))));
% 
