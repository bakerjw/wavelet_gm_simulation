% Program for generating artificial ground motions by stochastic ground
% motion model using wavelet packets
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

clear all;
clc;

% Rrup:     Rupture distance (km)
% Rhyp:     Hypocentral distance (km)
% inVs30:   Average shear wave velocity with 30m in surface
% Mw:       Minimum Moment Magnitude
% nsmpl:    # of samples to generate
% iflg:     =1  1 realization
%           =2  wave from median parameters
% offset:   # of offset for output filename

% data file that contains the result of regression analysis
filereg='.//regdata//regression_equation.txt';

iflg = 1;
offset=0;

% % conditions for each case
ind=1;
Mw(ind)=7;	Rrup(ind)= 10.0000;	Rhyp(ind)= 10.0000;	Vs30(ind)= 400.0;	nsmpl(ind)= 3;	ind=ind+1;
% Mw(ind)=8;	Rrup(ind)= 10.0000;	Rhyp(ind)= 10.0000;	Vs30(ind)= 400.0;	nsmpl(ind)= 3;	ind=ind+1;

% total of case 
ncase=ind-1;

% initializing random number generater
rs.TotalSign = RandStream('mt19937ar');
rs.MinorRand = RandStream('mt19937ar');
rs.MajorLoc  = RandStream('mt19937ar');
rs.MajorAmp  = RandStream('mt19937ar');
inrand=ceil(rand(1)*100);
for i=1:1:inrand;
    rand(rs.TotalSign);
    rand(rs.MinorRand);
    rand(rs.MajorLoc);
    rand(rs.MajorAmp);
end

% loop for each case
for j=1:1:ncase;
    % making directory for output files
    filename=sprintf('M%03.1f_Rr%08.4f_Rh%08.4f_Vs%06.1f',Mw(j),Rrup(j),Rhyp(j),Vs30(j));
    system(sprintf('mkdir %s',filename));

    % loop for each sample
    for i=1:1:nsmpl(j)
        disp(['computing... smpl' num2str(i,'% 6d') '/' num2str(nsmpl(j),'% 6d') '| case' num2str(j,'% 6d') '/' num2str(ncase,'% 6d')]);
        % generating each sample
        [th dt] = fn_get1Sim(Mw(j),Rrup(j),Vs30(j),rs,iflg,filereg,Rhyp(j));
        % storing each sample in NGA format
        fn_writeACCinNGA(sprintf('.\\%s\\%s_#%07d.txt',filename,filename,(i+offset)),length(th),dt,th, Mw(j), Rrup(j), Rhyp(j), Vs30(j))
    end;
    
end;

