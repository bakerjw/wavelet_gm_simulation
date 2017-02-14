function prm=fn_setParam
% setting control parameters
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

prm.general.gravity=980.665;

prm.general.compratio=0.70; % threshold of the ratio of major coefficients

prm.resp.hdamp=0.05;
% prm.resp.psuedo=1; % =0;psuedo vel from dis, =1;relative vel, =2;psuedo vel from acc
prm.resp.psuedo=0; % =0;psuedo vel from dis, =1;relative vel, =2;psuedo vel from acc

prm.epsiron.imethod=0;

prm.graph.flim=[0.3 50];
prm.graph.plim=[0.02 10];
prm.graph.flima=[0.0001 1];
prm.graph.plstime=0;
prm.graph.pletime=40;
prm.graph.xlimvaltrp=[0.01 10];
prm.graph.ylimvaltrp=[1 1000];

prm.general.unitfr = 0.2;    % 2*minimum frequency unit
prm.wavelet.wvlt   ='dmey';  % mother wavelet

totaldur = 10 * 2^14;   %163.84s = 0.01 * 2^14

prm.general.basedt         = 10; % prm.general.basedt = 0.010s
prm.general.ndata = round(totaldur/prm.general.basedt);
prm.wavelet.npow2 = round(log2(prm.general.ndata));


prm.general.basedt         = prm.general.basedt/1000;

%determine depth based on demand of frequency band

prm.general.mfold = prm.general.ndata/2 + 1;
prm.general.dur   = prm.general.basedt*prm.general.ndata;
prm.general.time  = prm.general.basedt:prm.general.basedt:prm.general.dur;
prm.general.nyqfr = round(1/(2.0*prm.general.basedt));
prm.general.df    = prm.general.nyqfr/prm.general.ndata * 2;
prm.general.freq  = 0:prm.general.df:prm.general.nyqfr;

dy=prm.general.nyqfr;
for i=15:-1:1
    dy=dy/2;
    if dy<prm.general.unitfr; depth=16-i; break; end;
end
prm.wavelet.depth=depth;                             % depth
prm.wavelet.dx=2^depth*prm.general.basedt;            % time resolution
prm.wavelet.dy=dy;                                   % frequency resolution
prm.wavelet.ord=fn_getORD(depth);                       % the order of the end nodes of wavelet packet tree
prm.wavelet.nTimeWindow=2^(prm.wavelet.npow2-depth);             % # of data in each end node
prm.wavelet.nFreqBand=2^depth;             % # of data in each end node
prm.wavelet.offset=12;             % offset of wavelet packet coefficients in time axis
