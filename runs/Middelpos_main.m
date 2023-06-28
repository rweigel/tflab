clear;

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

% Get input/output data
[B,E,t,infile,outfile] = Middelpos_clean(); 

if 1
    %% Band pass
    addpath(fullfile(scriptdir(),'..','fft'));
    Tm = 3*86400;
    band = [1/Tm,0.5];
    B = bandpass_(B,band);
    E = bandpass_(E,band);
    E = E(Tm+1:end-Tm,:);
    B = B(Tm+1:end-Tm,:);
end

%% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

if 1 
% Trim for faster run
B = B(1:10*pps,:);
E = E(1:10*pps,:);
t = t(1:10*pps,:);
end

%% Set output file base name using start/stop times of input data
filestr = sprintf('Middelpos-%s-%s',...
                  datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));

%% Set metadata
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = 1;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

%% Compute first TF
desc1 = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;

TF1 = tflab(B(:,1:2),E,opts1);
% Modify default description of run
TF1.Options.description = desc1;
TF1.Metadata = meta;
%TF1 = tflab_uncertainty(TF1);

fname1 = fullfile(scriptdir(),'data','Middelpos',[filestr,'-tf1.mat']);
savetf(TF1, fname1);

%% Compute second TF
desc2 = sprintf('OLS; One %d-day segment',size(B,1)/pps);
opts2 = tflab_options(1);
    opts2.tflab.loglevel = 1;

TF2 = tflab(B(:,1:2),E,opts2);
% Modify default description of run
TF2.Options.description = desc2;
TF2.Metadata = meta;

% Test S2.Z on same segments as S1.
%TF2 = tflab_metrics(TF2,opts2,TF1.Segment.IndexRange);
%TF2 = tflab_uncertainty(TF2);

fname2 = fullfile(scriptdir(),'data','Middelpos',[filestr,'-tf2.mat']);
savetf(TF2, fname2);

%Middelpos_plot;

%% Third TF
if (0)
    %%
    % Third TF
    desc3 = sprintf('OLS; %d 5-day segments',size(B,1)/(5*ppd));
    opts3 = tflab_options(1,dopts);
        opts3.tflab.loglevel = 1;
        opts3.td.window.width = 5*86400;
        opts3.td.window.shift = 5*86400;
        opts3.filestr = sprintf('%s-tf3',filestr);

    S3 = tflab(B(:,1:2),E,opts3);
    % Modify default description of run
    S3.Options.description = desc3;
    S3 = tflab_uncertainty(S3);

    fname = fullfile(scriptdir(),'data','Middelpos',[opts3.filestr,'.mat']);
    fprintf('Saving: %s\n',fname);
    save(fname,'-v7.3','-struct','S3');
    fprintf('Saved: %s\n',fname);
end