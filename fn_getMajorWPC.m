function [WpcMajorSim ierr]= fn_getMajorWPC(rs,xin,yin,prmcoef1,boundary)
% generating wavelet packets of major group
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% getMajorWPC returns wavelet packet coefficients in Major distribution
% WpcMajorSim   : wavelet packet coefficients in Major distribution
% ampMajorSim   : SQUARED amplitudes of wavelet packet coefficients in
% Major distirbution

minorElx =prmcoef1.minorElx;
minorEly =prmcoef1.minorEly;
minorSlx =prmcoef1.minorSlx;
minorSly =prmcoef1.minorSly;
minorRlxly=prmcoef1.minorRlxly;

majorEa  = prmcoef1.majorEa;
majorElx = prmcoef1.majorElx;
majorEly = prmcoef1.majorEly;
majorVlx = prmcoef1.majorVlx;
majorVly = prmcoef1.majorVly;
majorSlx = prmcoef1.majorSlx;
majorSly = prmcoef1.majorSly;
majorRlxly=prmcoef1.majorRlxly;

nxdim=length(xin);
nydim=length(yin);
ndata=nxdim*nydim;

y=(([0 yin(1:nydim-1)]+yin(1:nydim))./2);
x=(([0 xin(1:nxdim-1)]+xin(1:nxdim))./2);

% % generating amplitude of wavelet packets of major group
[ampMajorSim,nmajor]=fn_getMajorAmp(rs,majorEa,boundary,ndata);


% generating location using bivariate lognormal distribution
clear majorLLCov;
cov12=majorRlxly*majorSlx*majorSly;
majorLLCov=[majorVlx cov12; cov12 majorVly];

sfxyMajorSimt=zeros(nmajor,2);
xytmp=zeros(1,2);
ixy=0;
uxmx0=log(max(x));

uxmn0=log(min(x));
uymx0=log(max(y));
uymn0=log(min(y));

uxmx9=exp(uxmx0);
uxmn9=exp(uxmn0);
uymx9=exp(uymx0);
uymn9=exp(uymn0);

uxmx0=uxmx0-majorElx;
uxmn0=uxmn0-majorElx;
uymx0=uymx0-majorEly;
uymn0=uymn0-majorEly;

% % generating time and frequency location of wavelet packets of major group
while(true);
    ixy0=0;
    ixy=0;
    while(ixy0<nmajor);
        ixy=ixy+1;
        % using Box-Muller transform
        u1=rand(rs.MajorLoc);
        u2=rand(rs.MajorLoc);
        xytmp(1,1)=sqrt(-2.*log(u1)).*cos(2.*pi.*u2);
        xytmp(1,2)=sqrt(-2.*log(u1)).*sin(2.*pi.*u2);
        xytmp = (chol(majorLLCov,'lower')*xytmp')';
        ux = xytmp(1,1);
        uy = xytmp(1,2);

%         conditional mean and standard deviation on time axis for stopping
%         time
        mux_y=minorElx+minorRlxly*minorSlx*((majorEly + uy)-minorEly)/minorSly - majorElx;
        sgx_y=minorSlx*sqrt(1-minorRlxly^2);
        
        uxmx=min(uxmx0,mux_y+1*sgx_y);
        uxmn=uxmn0;
        uymx=uymx0;
        uymn=uymn0;
        
        if ixy>5000;    break;  end;
        
%         boundary for truncation
        if ux>uxmx; continue; end;
        if ux<uxmn; continue; end;
        if uy>uymx; continue; end;
        if uy<uymn; continue; end;

        ixy0=ixy0+1;
        sfxyMajorSimt(ixy0,1) = majorElx + ux;
        sfxyMajorSimt(ixy0,2) = majorEly + uy;
    end;
        clear WpcMajorSim;
        WpcMajorSim=zeros(nydim,nxdim);
    ierr=0;
    if ixy>5000;
        ierr=1;
        break;
    end;

        sfxMajorSim=zeros(1,nmajor);
        sfxMajorIndx=zeros(1,nmajor);
        sfyMajorSim=zeros(1,nmajor);
        sfyMajorIndx=zeros(1,nmajor);

        jdist=0;
        y2=y;
        y2(1)=y(1) + y(1)/2;

        nxdim2=round(nxdim/2)+1;
        x2=x;
        x2(1:(nxdim2-1))=0;
        x2(nxdim2:end)=((x([1:2:nxdim])+x([2:2:nxdim]))./2);
        
        for i=1:1:nmajor;
            tmpmdist=99999;
            for j=1:1:nydim;
                tmpdist=abs(exp(sfxyMajorSimt(i,2))-y2(j));
                if tmpdist<tmpmdist;
                    tmpmdist=tmpdist;
                    jdist=j;
                end;
            end;
            sfyMajorSim(i)=y2(jdist);
            sfyMajorIndx(i)=jdist;
            
            
            if jdist==1;
                xtt=x2;
                nxst=nxdim2;
            else
                xtt=x;
                nxst=1;
            end;
            
            jdist=0;

            tmpmdist=99999;
            for j=nxst:1:nxdim
                tmpdist=abs(exp(sfxyMajorSimt(i,1))-xtt(j));
                if tmpdist<tmpmdist;
                    tmpmdist=tmpdist;
                    jdist=j;
                end
            end
            sfxMajorSim(i)=xtt(jdist);
            sfxMajorIndx(i)=jdist;
            
            
        end;

        
        
        tamp = sum(ampMajorSim);
        z = ampMajorSim./tamp;
        ampMajorSim=z.*tamp;

        
        for i=1:1:nmajor
            if(WpcMajorSim(sfyMajorIndx(i),sfxMajorIndx(i))>10^(-5))
                WpcMajorSim(sfyMajorIndx(i),sfxMajorIndx(i))=WpcMajorSim(sfyMajorIndx(i),sfxMajorIndx(i))+ampMajorSim(i);
            else
                WpcMajorSim(sfyMajorIndx(i),sfxMajorIndx(i))=ampMajorSim(i);
            end
        end
        
    % common part
    WpcMajorSim=sqrt(WpcMajorSim);
    
if max(sfxMajorSim)>uxmx9; continue; end;
if min(sfxMajorSim)<uxmn9; continue; end;
if max(sfyMajorSim)>uymx9; continue; end;
if min(sfyMajorSim)<uymn9; continue; end;

    break;
    
end



