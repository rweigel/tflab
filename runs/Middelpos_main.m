% Get input/output data
[B,E,t,infile,outfile] = Middelpos_clean(); 

B = B(1:12*86400,:);
E = E(1:12*86400,:);
t = t(1:12*86400);

%% Band pass
Tm = 3*86400;
band = [1/Tm,0.5];
B = bandpass(B,band);
E = bandpass(E,band);
E = E(Tm+1:end-Tm,:);
B = B(Tm+1:end-Tm,:);

%%
% Make length an integer number of days.
ppd = 86400;
I = ppd*floor(size(B,1)/ppd);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

startstr = datestr(t(1),'yyyymmddTHHMMSS.FFF');
stopstr = datestr(t(end),'yyyymmddTHHMMSS.FFF');

filestr = sprintf('Middelpos-%s-%s',startstr,stopstr);

% Variable name information
iopts = struct('info',struct(),'td',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/km';
iopts.info.timeunit = 's';
iopts.info.timedelta = 1;
iopts.info.timestart = startstr;
iopts.info.stationid = 'Middelpos';
iopts.info.chainid = '';

%% Compute transfer functions

%%
% First TF
desc1 = sprintf('OLS; %d 1-day segments',size(B,1)/ppd);
opts1 = transferfnFD_options(1,iopts);
    opts1.transferfnFD.loglevel = 1;
    opts1.td.window.width = 86400;
    opts1.td.window.shift = 86400;
    opts1.filestr = sprintf('%s-tf1',filestr);

S1 = transferfnFD(B(:,1:2),E,opts1);
% Modify default description of run
S1.Options.description = desc1;
S1 = transferfnFDUncertainty(S1);

fname = fullfile(scriptdir(),'data','Middelpos',[opts1.filestr,'.mat']);
fprintf('Saving: %s\n',fname);
save(fname,'-v7.3','-struct','S1');
fprintf('Saved: %s\n',fname);

%%
% Second TF
desc2 = sprintf('OLS; One %d-day segment',size(B,1)/ppd);
opts2 = transferfnFD_options(1,iopts);
    opts2.transferfnFD.loglevel = 1;
    opts2.filestr = sprintf('%s-tf2',filestr);

S2 = transferfnFD(B(:,1:2),E,opts2);           
% Modify default description of run
S2.Options.description = desc2;

% Test S2.Z on same segments as S1.
S2 = transferfnFDMetrics(S2,opts2,S1.Segment.IndexRange);
S2 = transferfnFDUncertainty(S2);

fname = fullfile(scriptdir(),'data','Middelpos',[opts2.filestr,'.mat']);
fprintf('Saving: %s\n',fname);
save(fname,'-v7.3','-struct','S2');
fprintf('Saved: %s\n',fname);

%%
% Third TF
desc3 = sprintf('OLS; %d 5-day segments',size(B,1)/(5*ppd));
opts3 = transferfnFD_options(1,iopts);
    opts3.transferfnFD.loglevel = 1;
    opts3.td.window.width = 5*86400;
    opts3.td.window.shift = 5*86400;
    opts3.filestr = sprintf('%s-tf3',filestr);

S3 = transferfnFD(B(:,1:2),E,opts3);
% Modify default description of run
S3.Options.description = desc3;
S3 = transferfnFDUncertainty(S3);

fname = fullfile(scriptdir(),'data','Middelpos',[opts3.filestr,'.mat']);
fprintf('Saving: %s\n',fname);
save(fname,'-v7.3','-struct','S3');
fprintf('Saved: %s\n',fname);
