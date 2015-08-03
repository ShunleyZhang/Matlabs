I = load('cost_0.txt');
I2 = load('cost_1.txt');
P = load('PatchVector14.txt');
D = load('pair_0.txt');

%%
% [ GPMap ] = GoodPixelSelect( I, I2, 0.3 );
% LRvec = zeros(size(GPMap,1), size(P,2)+1);
% idx = 1;
% tmp = 1;
% s = size(P,2);
% for i=1:size(I,1)
%   for j = 1:size(I,2)
%       if GPMap(i,j) == 1
%           for a = 1 : s
%              LRvec(tmp,a) = P(idx,a);
%           end
%           LRvec(tmp, s+1) = D(i,j);
%       tmp = tmp + 1;
%       end
%       idx = idx + 1;
%   end 
% end
% 
% 
% fid = fopen('GPMap.txt', 'w');
% for i=1:size(I,1)
%   for j = 1:size(I,2)
%       fprintf(fid, '%d ', GPMap(i,j));
%   end 
%   fprintf(fid, '\n');
% end
% fclose(fid);
% 
% fid = fopen('LRvec.txt', 'w');
% for i=1:size(LRvec,1)
%   for j = 1:size(LRvec,2)
%       fprintf(fid, '%d ', LRvec(i,j));
%   end 
%   fprintf(fid, '\n');
% end

LRvec = load('LRvec.txt');
disp(size(LRvec));
[A_hat,E_hat] = exact_alm_rpca(LRvec, 0.01);
fid = fopen('LRvecR.txt', 'w');
for i=1:size(LRvec,1)
  for j = 1:size(LRvec,2)
      fprintf(fid, '%.0f ', A_hat(i,j));
  end 
  fprintf(fid, '\n');
end

A_hat = load('LRvecR.txt');
E = D;

tmp = 1;
for i=1:size(I,1)
  for j = 1:size(I,2)
      if GPMap(i,j) == 1
          D(i,j) = A_hat(tmp, s+1);
          
          su = E_hat(tmp,:);
            disp(su);
          if sum(su(:)) > 5
              disp(sum(su));
              D(i,j) = 0;
          end
          tmp = tmp+1;
      else
          D(i,j) = 0;
          E(i,j) = 0;
      end
  end 
end
X = D(:,26) - E(:,26);
disp(max(X));
x = 1:size(I,1);
y = 1:size(I,2);
[X,Y] = meshgrid(x,y);

%surf(X,Y,I');
subplot(211);
imagesc(E);
subplot(212);
imagesc(D);