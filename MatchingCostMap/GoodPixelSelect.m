function [ GPMap ] = GoodPixelSelect( Cost1, Cost2, DeltaTh )
%GOODPIXELSELECT 此处显示有关此函数的摘要
%   此处显示详细说明
%   Cost1 Cost2:    Minimum and second minimum Matching Cost
%   DeltaTh:        Delta Threshold
%   PixelList:      Map of Good Pixel, 1 for good pixel, 0 for the others

s = size(Cost1);
GPMap = zeros(s(1), s(2));
for i=1:s(1)
  for j = 1:s(2)
      if Cost2(i,j) > Cost1(i,j)*(1+DeltaTh)
          GPMap(i,j)=1;
      end
  end 
end
end

