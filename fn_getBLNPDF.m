function outPDF=fn_getBLNPDF(xin,yin,Elx,Ely,Vlx,Vly,lxlyMinorCor)
% generating probability density function of bivariate lognormal distribution
% truncated by stoppin time for each frequency band
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

nxdim=length(xin);
nydim=length(yin);

% y=exp((log(yin(1:nydim))+log([df yin(1:(nydim-1))]))./2);
y=(([0 yin(1:nydim-1)]+yin(1:nydim))./2);
x=(([0 xin(1:nxdim-1)]+xin(1:nxdim))./2);

Slx=sqrt(Vlx);
Sly=sqrt(Vly);

outPDF=zeros(nydim,nxdim);
rho2= 1 - lxlyMinorCor^2;

nxdim2=round(nxdim/2) + 1;
x2=x;
x2(1:(nxdim2-1))=0;
x2(nxdim2:end)=((x([1:2:nxdim])+x([2:2:nxdim]))./2);
y2=y;
y2(1)=y(1) + y(1)/2;
x_prime = log(x) - Elx;
x_prime2 = log(x2) - Elx;
y_prime2 = log(y2) - Ely;

for i=1:1:nydim;

    nxdim3=1;
    x3=x;
    x_prime3=x_prime;
    if i==1;
        nxdim3=nxdim2;
        x3=x2;
        x_prime3=x_prime2;
    end;

    mux_y=Elx+lxlyMinorCor*Slx*(y2(i)-Ely)/Sly;
    sgx_y=Slx*sqrt(1-lxlyMinorCor^2);
    bnd=exp(mux_y+2*sgx_y);

    for j=nxdim3:1:nxdim;
        if j>nxdim;
            continue;
        end;
        if bnd<x3(j)
            continue;
        end;

        outPDF(i,j)=  1/( 2*pi*Slx*Sly*sqrt(rho2) ) /x3(j) /y2(i)...
     *exp( -1/(2*rho2) * (x_prime3(j)^2./Vlx + y_prime2(i)^2./Vly - 2*lxlyMinorCor*x_prime3(j)*y_prime2(i)/Slx/Sly) );
    end;
end;

outPDF=outPDF./sum(sum(outPDF(:,:)));


