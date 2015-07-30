I = load('Depth_comp_rstd.txt');
s = size(I);
fp = fopen('Depth_compaV_03.txt','wt');

    for j = 1:s(1)
        for k = 1:s(2)
                I(j,k) = min(I(j,k)*3, 1);
                fprintf(fp,'%.0f ',I(j,k));
         end 
        fprintf(fp , '\n');
    end