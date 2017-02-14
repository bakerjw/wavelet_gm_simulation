function ord=fn_getORD(depth)
% getting the order of wavelet packet nodes in frequency
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

ndim=2^depth;
ord=[1:1:ndim];
for i=1:1:depth-1
    ordt=ord;
    mdim1=ndim/(2^i);
    mdim2=mdim1/2;
    mdim0=mdim1;
    for j=1:1:2^(i-1)
        ord(mdim0+1:mdim0+mdim2)=ordt(mdim0+mdim2+1:mdim0+2*mdim2);
        ord(mdim0+mdim2+1:mdim0+2*mdim2)=ordt(mdim0+1:mdim0+mdim2);
        mdim0=mdim0+2*mdim1;
    end
end

