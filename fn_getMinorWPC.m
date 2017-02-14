function [WpcMinorSimR]= fn_getMinorWPC(rs,x,y,prmcoef1,wpcMinorEnergy)
% generating wavelet packets of minor group
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% generating probability density function of bivariate lognormal
% distribution
WpcMinorSimNR0 = fn_getBLNPDF(x,y, ...
    prmcoef1.minorElx,prmcoef1.minorEly,prmcoef1.minorVlx,prmcoef1.minorVly,prmcoef1.minorRlxly);

nxdim=size(WpcMinorSimNR0,2);
nydim=size(WpcMinorSimNR0,1);

% random factor of wavelet coefficients of minor group
aa=randn(rs.MinorRand,nydim,nxdim);
WpcMinorSimR= exp(aa.*prmcoef1.minorRnd).* WpcMinorSimNR0;

% normalize to conserve energy
WpcMinorSimR=sqrt(wpcMinorEnergy.*WpcMinorSimR./sum(sum(WpcMinorSimR)));

