function [ampMajorSim,nmajor]=fn_getMajorAmp(rs,majorEa,boundary,ndata)
% generating amplitude of wavelet packets in major group
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

iamp=0;
iamp0=0;
tmpsum=0;
ampMajorSimTmp=zeros(1,ndata);

while(tmpsum<boundary);
    iamp=iamp+1;
%         inverse methog from uniform to exponential
    u1 = -log(rand(rs.MajorAmp,1,1))*majorEa;

    iamp0=iamp0+1;
    ampMajorSimTmp(iamp0) = u1;

    tmpsum=tmpsum+ampMajorSimTmp(iamp0);
%     going out when energy goes beyond its bounday
end;

nmajor=iamp0;
if(abs(tmpsum-boundary)>abs(tmpsum-ampMajorSimTmp(iamp)-boundary))
    nmajor=nmajor-1;
end
ampMajorSim =zeros(1,nmajor);
for i=1:1:nmajor
    ampMajorSim(i)=ampMajorSimTmp(i);
end
clear ampMajorSimTmp;

