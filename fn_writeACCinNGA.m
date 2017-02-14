function fn_writeACCinNGA(filename,ndata,dt,acc, M, Rrup, Rhyp, Vs30)
% storing time history data in NGA format
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

fid = fopen(filename,'w');
fprintf(fid,'%s\n','Stochastic Ground Motion By Yoshi');
fprintf(fid,'M=%-10.4f / Rrup=%-10.2fkm / Rhyp=%-10.2fkm / Vs30=%-10.3fm/s\n',M,Rrup,Rhyp,Vs30);
fprintf(fid,'%s\n','ACCELERATION TIME HISTORY IN UNITS OF G');
fprintf(fid,'%-10d%-10.6f%s\n',ndata,dt,'    NPTS, DT');
fprintf(fid,'%-15.6E%-15.6E%-15.6E%-15.6E%-15.6E\n',acc);
fclose(fid);

