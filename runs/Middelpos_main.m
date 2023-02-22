clear;

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

% Get input/output data
[B,E,t,infile,outfile] = Middelpos_clean(); 

if 1 % For testing a run that takes less time.
    B = B(1:12*86400,:);
    E = E(1:12*86400,:);
    t = t(1:12*86400);
end

if 1
%% Band pass
    addpath(fullfile(scriptdir(),'..','fft'));
    Tm = 3*86400;
    band = [1/Tm,0.5];
    B = bandpass(B,band);
    E = bandpass(E,band);
    E = E(Tm+1:end-Tm,:);
    B = B(Tm+1:end-Tm,:);
end

%% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

%% Set output file base name using start/stop times of input data
filestr = sprintf('Middelpos-%s-%s',...
                datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));


%% Set variable name information
iopts = struct('info',struct(),'td',struct());
iopts.info.instr     = {'$B_x$','$B_y$'};
iopts.info.inunit    = 'nT';
iopts.info.outstr    = {'$E_x$','$E_y$'};
iopts.info.outunit   = 'mV/km';
iopts.info.timeunit  = 's';
iopts.info.timedelta = 1;
iopts.info.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
iopts.info.chainid   = ''; % Use SANSA?
iopts.info.stationid = 'Middelpos';


%% Compute first TF
desc1 = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts1 = transferfnFD_options(1,iopts);
    opts1.transferfnFD.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;
    opts1.filestr = sprintf('%s-tf1',filestr);

S1 = transferfnFD(B(:,1:2),E,opts1);
% Modify default description of run
S1.Options.description = desc1;
S1 = transferfnFD_uncertainty(S1);

fname1 = fullfile(scriptdir(),'data','Middelpos',[opts1.filestr,'.mat']);
savetf(S1, fname1);


%% Compute second TF
desc2 = sprintf('OLS; One %d-day segment',size(B,1)/pps);
opts2 = transferfnFD_options(1,iopts);
    opts2.transferfnFD.loglevel = 1;
    opts2.filestr = sprintf('%s-tf2',filestr);

S2 = transferfnFD(B(:,1:2),E,opts2);           
% Modify default description of run
S2.Options.description = desc2;

% Test S2.Z on same segments as S1.
S2 = transferfnFD_metrics(S2,opts2,S1.Segment.IndexRange);
S2 = transferfnFD_uncertainty(S2);

fname = fullfile(scriptdir(),'data','Middelpos',[opts2.filestr,'.mat']);
savetf(S2, fname);


%% Third TF
if (0)
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
    S3 = transferfnFD_uncertainty(S3);

    fname = fullfile(scriptdir(),'data','Middelpos',[opts3.filestr,'.mat']);
    fprintf('Saving: %s\n',fname);
    save(fname,'-v7.3','-struct','S3');
    fprintf('Saved: %s\n',fname);
end