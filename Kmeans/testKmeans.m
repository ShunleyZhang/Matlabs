y = load('PatchVector14.txt');
disp(size(y));
I = zeros(271, 404); %414 281    429 378   574 397
%调用Kmeans函数
%X N*P的数据矩阵
%Idx N*1的向量,存储的是每个点的聚类标号
%Ctrs K*P的矩阵,存储的是K个聚类质心位置
%SumD 1*K的和向量,存储的是类间所有点与该类质心点距离之和
%D N*K的矩阵，存储的是每个点与所有质心的距离;
 
%Idx,Ctrs,SumD,D] = kmeans(X,3,'Replicates',3,'Options',opts);
[Idx,Ctrs,SumD,D] = kmeans(y,8);
A = [0,0,0,0,0,0,0,0];
A2 = zeros(1,16);
for i = 1:size(y,1)
    tmp = Idx(i);
    A(tmp) = A(tmp)+1;
end
disp(A);
%find A max
[Amax,maxId] = max(A);
disp(Amax);
disp(maxId);
B = zeros(Amax,1);
s = size(y);
C = zeros(Amax,s(2));
tmp = 1;
id = 1;
% fid = fopen('cluster.txt','w');
for i = 1:s(1)
    if(Idx(id) == maxId)
        C(tmp,:) = y(id,:);
        tmp = tmp + 1;
%         for j = 1:size(y,2)
%             fprintf(fid, '%d ', y(id,j));
%         end
%         fprintf(fid,'\n');
    end
    id = id + 1;
end
disp(size(C));
[Idx2,Ctrs,SumD,D] = kmeans(C,16);
for i = 1:size(C,1)
    tmp = Idx2(i);
    A2(tmp) = A2(tmp)+1;
end
disp(A2);
%find A max
[Amax,maxId] = max(A2);
disp(Amax);
disp(maxId);
tmp = 1;
id = 1;
for i = 1:size(y,1)
    if(Idx(id) == maxId)
        Idx(id) = Idx2(tmp)/4 + 8;
        tmp = tmp + 1;
    end
    id = id+ 1;
end
index = 1;
for i = 1:271
    for j = 1:404
        I(i,j) = Idx(index);
       % disp(I(i,j));
        index = index+1;
    end
end

x = 1:size(I,1);
y = 1:size(I,2);
[X,Y] = meshgrid(x,y);

%surf(X,Y,I');

imagesc(I); colormap([1 1 1;1 1 0;1 0 0;1 0 1; 0 0 1; 0 1 1; 0 1 0; 0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2]);colorbar;
% %409 276
% % [Idx,Ctrs,SumD,D] = kmeans(y,4);
% % 
% % disp(Ctrs);
% % 
% % bandwidth = 0.5;
% % tic
% % [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(y',bandwidth);
% % toc
% % disp(clustCent);


 
