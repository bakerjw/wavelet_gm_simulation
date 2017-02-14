function [outprm]=fn_getPARAM8(Mw,Rrup,Rhyp,Vs30,filereg)
% generating parameters based on the regression equations
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu


if (nargin < 5)
 disp('need more input');
 stop;
end

D=Rrup;
% initialize outprm
outprm = fn_initDBs(2);

% generate parameter with corrrelation of mulivariate normal distribution

% Mw: Moment Magnitude
% Dhyp: Hypocentral Distance(km)
% Drup: Closest Distance(km)
% Vs30: Average of shear wave velocity(m/s)
% filereg: name of file that contains regression equations

lD=log(D);
lVs30=log(Vs30);

filename=char(filereg);

% reading regression equations
    [acoef,bcoef,ccoef,dcoef,ecoef,fcoef,epsEi,epsSi, ...
        covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9,covi10,covi11,covi12, ...
        epsEm,epsSm, ...
        covm1,covm2,covm3,covm4,covm5,covm6,covm7,covm8,covm9,covm10,covm11,covm12,mn13arr,sn13arr ...
        ] = ...
        textread(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ,'delimiter', '\t', 'emptyvalue', NaN,'headerlines',1,'whitespace', '','endofline','\r\n');
    covi=[covi1 covi2 covi3 covi4 covi5 covi6 covi7 covi8 covi9 covi10 covi11 covi12];
    covm=[covm1 covm2 covm3 covm4 covm5 covm6 covm7 covm8 covm9 covm10 covm11 covm12];

    D2=Rhyp-Rrup;
    D1=D;
    lAVS30=lVs30;
    
    mn13=mn13arr(1);
    sn13=sn13arr(1);
    m13=log(mn13^2/sqrt(sn13+mn13^2));
    s13=sqrt(log(sn13/mn13^2+1));

%     median prediction of each parameter
    ind=1;
    hcoef=10;
    hcoef1=1;
    majorExM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*exp(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef1^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    majorEyM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    majorSxM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*exp(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef1^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    majorSyM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    majorRxyM = acoef(ind) + bcoef(ind)*Mw + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;

    minorExM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*exp(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef1^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    minorEyM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    minorSxM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*exp(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef1^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    minorSyM  = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    minorRxyM = acoef(ind) + bcoef(ind)*Mw + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    
    majorEaM     = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;
    totalEnergyM = acoef(ind) + bcoef(ind)*Mw + ccoef(ind)*log(Mw) + dcoef(ind)*D2 + ecoef(ind)*log(sqrt(D1^2+hcoef^2)) + fcoef(ind)*lAVS30; ind=ind+1;

        while true
%             intra-event residuals
            epsi = mvnrnd(epsEi,covi);
%             inter-event residuals
            epsm=mvnrnd(epsEm,covm);
%             total residuals
            epst=epsi+epsm;
%             random factor of wavelet packets of minor group
            minorRnd=lognrnd(m13,s13);
            
%     prediction of each parameter with residuals
            ind=1;
            majorEx = exp(majorExM + epst(ind)); ind=ind+1;
            majorEy = exp(majorEyM + epst(ind)); ind=ind+1;
            majorSx = exp(majorSxM + epst(ind)); ind=ind+1;
            majorSy = exp(majorSyM + epst(ind)); ind=ind+1;
            majorRxy = 2*normcdf(majorRxyM + epst(ind),0,1)-1; ind=ind+1;

            minorEx = exp(minorExM + epst(ind)); ind=ind+1;
            minorEy = exp(minorEyM + epst(ind)); ind=ind+1;
            minorSx = exp(minorSxM + epst(ind)); ind=ind+1;
            minorSy = exp(minorSyM + epst(ind)); ind=ind+1;
            minorRxy = 2*normcdf(minorRxyM + epst(ind),0,1)-1; ind=ind+1;
            
            majorEa = exp(majorEaM + epst(ind)); ind=ind+1;
            totalEnergy = exp(totalEnergyM + epst(ind)); ind=ind+1;
            
            majorVx = majorSx^2;
            majorVy = majorSy^2;
            majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;

            majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
            majorVlx = log(majorVx/majorEx^2 + 1);
            majorSlx = sqrt(majorVlx);
            majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
            majorVly = log(majorVy/majorEy^2 + 1);
            majorSly = sqrt(majorVly);

            majorCovlxly = log(majorExy/majorEx/majorEy);
            majorRlxly = majorCovlxly/majorSlx/majorSly;

            minorVx = minorSx^2;
            minorVy = minorSy^2;
            minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
            minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
            minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
            minorVlx = minorSlx^2;
            minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
            minorSly = sqrt(log(minorVy/minorEy^2 + 1));
            minorVly = minorSly^2;

            minorCovlxly = log(minorExy/minorEx/minorEy);
            minorRlxly = minorCovlxly/minorSlx/minorSly;
            
            majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
            [T,err] = cholcov(majorLLCov);
            if err == 0;
                break;
            end;
        end;

%         output parameters
        outprm(1) = struct('M',Mw,'hdist',D,'vs30',Vs30, ...
                'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
                'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
                'totalEnergy',totalEnergy, 'majorEa',majorEa,'minorRnd',minorRnd ...
            );

%         for median prediction
    majorExM = exp(majorExM);
    majorSxM = exp(majorSxM);
    majorEyM = exp(majorEyM);
    majorSyM = exp(majorSyM);
    majorRxyM = 2*normcdf(majorRxyM,0,1)-1;

    minorExM = exp(minorExM);
    minorSxM = exp(minorSxM);
    minorEyM = exp(minorEyM);
    minorSyM = exp(minorSyM);
    minorRxyM = 2*normcdf(minorRxyM,0,1)-1;

    majorEaM = exp(majorEaM);
    totalEnergyM = exp(totalEnergyM);

%     about major
    majorVxM = majorSxM^2;
    majorVyM = majorSyM^2;
    majorExyM = majorRxyM*majorSxM*majorSyM + majorExM*majorEyM;
    
    majorElxM = log(majorExM^2 / sqrt(majorVxM+majorExM^2));
    majorSlxM = sqrt(log(majorSxM/majorExM^2 + 1));
    majorVlxM = majorSlxM^2;
    majorElyM = log(majorEyM^2 / sqrt(majorVyM+majorEyM^2));
    majorSlyM = sqrt(log(majorSyM/majorEyM^2 + 1));
    majorVlyM = majorSlyM^2;
    
    majorCovlxlyM = log(majorExyM/majorExM/majorEyM);
    majorRlxlyM = majorCovlxlyM/majorSlxM/majorSlyM;
    
%     about minor
    minorVxM = minorSxM^2;
    minorVyM = minorSyM^2;
    minorExyM = minorRxyM*minorSxM*minorSyM + minorExM*minorEyM;

    minorElxM = log(minorExM^2 / sqrt(minorVxM+minorExM^2));
    minorSlxM = sqrt(log(minorSxM/minorExM^2 + 1));
    minorVlxM = minorSlxM^2;
    minorElyM = log(minorEyM^2 / sqrt(minorVyM+minorEyM^2));
    minorSlyM = sqrt(log(minorSyM/minorEyM^2 + 1));
    minorVlyM = minorSlyM^2;
    
    minorCovlxlyM = log(minorExyM/minorExM/minorEyM);
    minorRlxlyM = minorCovlxlyM/minorSlxM/minorSlyM;
    
        outprm(2) = struct('M',Mw,'hdist',D,'vs30',Vs30, ...
                'minorElx',minorElxM, 'minorSlx',minorSlxM, 'minorVlx',minorVlxM, 'minorEly',minorElyM, 'minorSly',minorSlyM, 'minorVly',minorVlyM,'minorRlxly',minorRlxlyM, ...
                'majorElx',majorElxM, 'majorSlx',majorSlxM, 'majorVlx',majorVlxM, 'majorEly',majorElyM, 'majorSly',majorSlyM, 'majorVly',majorVlyM,'majorRlxly',majorRlxlyM, ...
                'totalEnergy',totalEnergyM, 'majorEa',majorEaM, 'minorRnd',minorRnd ...
            );

    
