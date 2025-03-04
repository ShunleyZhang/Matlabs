Map = cell(4,4);
%%depth map
Map{1,1} = load('pair_0 (2).txt');%0
Map{1,2} = load('pair_0 (3).txt');%0.25
Map{1,3} = load('pair_0.txt');%0.5
Map{1,4} = load('pair_0 (4).txt');%0.75
%%cost1
Map{2,1} = load('cost_0 (2).txt');%0
Map{2,2} = load('cost_0 (3).txt');%0.25
Map{2,3} = load('cost_0.txt');%0.5
Map{2,4} = load('cost_0 (4).txt');%0.75
%%cost2
Map{3,1} = load('cost_1 (2).txt');%0
Map{3,2} = load('cost_1 (3).txt');%0.25
Map{3,3} = load('cost_1.txt');%0.5
Map{3,4} = load('cost_1 (4).txt');%0.75
for i = 1:4
    Map{3,i} = max(Map{3,i},0);
    Map{2,i} = max(Map{2,i},0);
    [Map{4,i}] = GoodPixelSelect( Map{2,i}, Map{3,i}, 0.3 );
end

for i = 1:4
    subplot(4,4,(i-1)*4+1);   imagesc(Map{1,i});
    subplot(4,4,(i-1)*4+2);   imagesc(Map{2,i});
    subplot(4,4,(i-1)*4+3);   imagesc(Map{3,i}./Map{2,i});
    subplot(4,4,(i-1)*4+4);   imagesc(Map{1,i}.*Map{4,i});
end
% x = 1:size(I,1);
% y = 1:size(I,2);
% [X,Y] = meshgrid(x,y);

%surf(X,Y,I');
%depth
% subplot(441);   imagesc(O);
% subplot(442);   imagesc(P);
% subplot(443);   imagesc(N);
% subplot(444);   imagesc(Q);
% %cost
% subplot(447);   imagesc(I);
% subplot(248);   imagesc(I2./I);
%cost2/cost1
%good pixel