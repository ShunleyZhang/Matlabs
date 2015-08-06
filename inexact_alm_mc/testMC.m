I = imread('wall.png');
I1 = im2double(I);
D = I1(:,:,1);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if D(i,j)>0.80
%             disp(D(i,j))
            D(i,j) = 0;
        end
    end
end

Ds = sparse(D);

[A,S] = inexact_alm_mc(Ds, 0.5);
U = diag(A.U);
V = diag(A.V);
S1 = zeros(size(U,2), size(V,1));
S1(1:size(S,1),1:size(S,1)) = S;

subplot(311);
imshow(full(A.U* A.V'));
% [A2] = inexact_alm_mc(Ds, 0.9);
% subplot(312);
% disp(svp);
% imshow(full(A2.U * svp * A2.V'));
% [A3] = inexact_alm_mc(Ds, 2);
% subplot(313);
% disp(svp);
% imshow(full(A3.U * svp * A3.V'));