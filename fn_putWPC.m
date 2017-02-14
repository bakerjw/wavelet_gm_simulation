function y=fn_putWPC(ord,wpcmat,depth,wvlt,dt)
% inverse wavelet packet transform
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% putWPC returns time history data from wavelet packet coefficients.
% ord       : order of wavelet packet nodes in wavelet tree
% wpcmat    : wavelet packet coefficients
% wvlt      : type of mother wavelet
% ndata     : # of data of time history

wpcmat=wpcmat./sqrt(dt);

sizeWpc = size(wpcmat);
ndata=sizeWpc(1)*sizeWpc(2);
nydim = sizeWpc(1);
nxdim = max(read(wpdec(zeros(1,ndata),depth,wvlt),'sizes',2^depth-1));

nTimeWindow=round(ndata*2^(-depth));             % # of data in each end node
tail=(nxdim-nTimeWindow)/2;
first= 1+floor(tail);
last = nxdim-ceil(tail);

wpcmattmp=zeros(1,nxdim*nydim);

tord=zeros(nydim,1);
for i=1:1:nydim
    stp=(ord(i)-1)*nxdim;
    wpcmattmp(1,stp+first:stp+last)=wpcmat(i,:);
    tord(i,1)=nydim+i-2;
end
twpc=cfs2wpt(wvlt,ndata,tord,2,wpcmattmp);


nlen=last-first+1;
nlen2=nlen/2;
twpc=wpsplt(twpc,[depth 0]);
tord1=max(tord)+1;
tord2=tord1+1;

nxdim2=max(read(twpc,'sizes',(2^(depth+1))));
effnx=round(ndata/(2^(depth+1)));
tail1=(nxdim2-effnx)/2;
first1= 1+floor(tail1);
last1 = nxdim2-ceil(tail1);

cfs1=zeros(1,nxdim2);
cfs2=zeros(1,nxdim2);
cfs2(first1:last1)=wpcmat(1,nlen2+1:nlen);

twpc=write(twpc,'cfs',tord1,cfs1,'cfs',tord2,cfs2);

y=wprec(twpc);