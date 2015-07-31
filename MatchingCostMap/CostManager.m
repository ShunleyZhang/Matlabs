I = load('cost_0.txt');
% I2 = load('cost_1.txt');
% M = load('pair_0.txt');
% N = load('pair_1.txt');
% I = max(I,0);I2=max(I2,0);
% 
% for i=1:size(I,1)
%   for j = 1:size(I,2)
%       if I2(i,j)<I(i,j)*1.5 
%           M(i,j)=0;
%       end
%   end 
% end

x = 1:size(I,1);
y = 1:size(I,2);
[X,Y] = meshgrid(x,y);

%surf(X,Y,I');
imagesc(I);