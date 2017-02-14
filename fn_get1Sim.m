function [th dt] = fn_get1Sim(Mw,Rrup,Vs30,rs,iflg,filereg,Rhyp)
% generating ONE ground motion
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% Mw      : Moment Magnitude
% Rhyp    : Hypocentral Distance (km)
% Rrup    : Closest Distance (km)
% Vs30    : Average shear wave velocity with 30m in surface
% rs      : handles of Random variables
% filereg : name of file that contains regression equations
% iflg    : =1  1 realization   : =2  wave from median parameters

if (nargin < 7)
%     if hypocentral distance is missing, hypocentral distance equals to
%     closest distance
    Rhyp = Rrup;
end

% set parameters
prm=fn_setParam();
ierr=1;
while(ierr~=0);
% get parameters from regression equations
    [prmcoef0]=fn_getPARAM8(Mw,Rrup,Rhyp,Vs30,filereg);

    prmcoef=prmcoef0(iflg);
    x=prm.wavelet.dx.*[1:1:prm.wavelet.nTimeWindow];
    y=prm.wavelet.dy.*[1:1:prm.wavelet.nFreqBand];

    totalEnergy = prmcoef.totalEnergy;
    minorEnergy=(1-prm.general.compratio)*totalEnergy;
    majorEnergy=prm.general.compratio*totalEnergy;

    nx=length(x);
    ny=length(y);

    % these values are from the computation of correlation of sign of the data
    % with available period contens in g.e. 10 sec  and  dt of wavelet packet
    % coefficients is 2.56 sec.
    yrand = randn(rs.TotalSign,ny,nx);

    % generate minor distribution
    [WpcMinorSimR] = fn_getMinorWPC(rs,x,y,prmcoef,minorEnergy);

    % generate major distribution
    [WpcMajorSim ierr] = fn_getMajorWPC(rs,x,y,prmcoef,majorEnergy);
    if ierr~=0;
        disp('recompute parameters again');
    end;
    
end;
% combining major and minor group
[WpcSim]=sqrt(WpcMinorSimR.*WpcMinorSimR+WpcMajorSim.*WpcMajorSim);

% one more decompose and set zero into the wavelet packets in the lowest
% frequency band
wpc2=zeros(size(WpcSim,1),size(WpcSim,2));
offset=prm.wavelet.offset;
wpc2(2:end,(offset+1):end) = WpcSim(2:end,1:(end-offset));
offset2=round(offset/2);
wpc2(1,(offset2+1):end) = WpcSim(1,1:(end-offset2));
WpcSim = wpc2;

% adding random sign
WpcSim(:,:)= sign(yrand(:,:)).*abs(WpcSim(:,:));

dt = prm.general.basedt;

% inverse wavelet packet transform
th = fn_putWPC(prm.wavelet.ord,WpcSim,prm.wavelet.depth,prm.wavelet.wvlt,dt);
mxdur=exp(prmcoef.minorElx+2*prmcoef.minorSlx)+offset*prm.wavelet.dx+30;
mxstp =min(ceil(mxdur/dt),prm.general.ndata);

th = th(1:mxstp);