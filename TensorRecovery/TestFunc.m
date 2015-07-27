% alpha = [1,1,1];
% beta = alpha;
% gamma = alpha/2;
% disp(gamma);

%%test SVD and UD(A)V
Mat = [
1 2 3
4 5 6
7 8 9
];
Mat(:,:)=0;
disp(Mat);
% % disp(S); disp(U);
% [U,S,V]=svd(Mat);
% diagS = diag(S);
% svp = length(find(diagS > 1/mu));
%  disp(svp);
%  disp(diagS);
% mu=1;
% A = diag(diagS(1:svp)-1/mu);
% afterD = max(diagS(1:svp)-1/mu, 0);
% afterD2 = max(Mat-1/mu, 0);
% disp(afterD2);
% disp(A);
%  disp(afterD);
% D3 = zeros(2,4,3);
% for i=1:2
% for j=1:4
% for k=1:3
% D3(i,j,k)=i+j+k;
% end
% end
% end
% dim = size(D3);
% i=3;
% X = reshape(shiftdim(D3,i-1), dim(i), []);
% disp(size(X));
% dim = circshift(dim, [1-i, 1-i]);
% X = shiftdim(reshape(X, dim), length(dim)+1-i);
% disp(size(X));