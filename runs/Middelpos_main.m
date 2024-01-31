clear;

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

% Get input/output data
[Bo,Eo,to,infile,outfile] = Middelpos_clean(); 

%% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(Bo,1)/pps);
Bo = Bo(1:I,:);
Eo = Eo(1:I,:);
to = to(1:I);

if 1
    % Trim for faster run
    Bo = Bo(1:10*pps,:);
    Eo = Eo(1:10*pps,:);
    to = to(1:10*pps,:);
end

%% Set output file base name using start/stop times of input data
filestr = sprintf('Middelpos-%s-%s',...
                  datestr(to(1),'yyyymmdd'),datestr(to(end),'yyyymmdd'));

%% Set metadata
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = 1;
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1;
    meta.timestart = datestr(to(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

if 1  
    %% First TF
    
    Tm = 6*86400;
    band = [1/Tm,0.5];
    B = bandpass_(Bo,band);
    E = bandpass_(Eo,band);
    %E = E(Tm+1:end-Tm,:);
    %B = B(Tm+1:end-Tm,:);

    desc1 = sprintf('6-day bandpass');
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
end

if 1
    %% Second TF

    %% Band pass
    Tm = 3*86400;
    band = [1/Tm,0.5];
    B = bandpass_(Bo,band);
    E = bandpass_(Eo,band);
    %E = E(Tm+1:end-Tm,:);
    %B = B(Tm+1:end-Tm,:);

    desc2 = sprintf('3-day bandpass');
    opts2 = tflab_options(1);
        opts2.tflab.loglevel = 1;
        opts2.td.window.width = pps;
        opts2.td.window.shift = pps;

    TF2 = tflab(B(:,1:2),E,opts2);
    % Modify default description of run
    TF2.Options.description = desc2;
    TF2.Metadata = meta;

    fname1 = fullfile(scriptdir(),'data','Middelpos',[filestr,'-tf2.mat']);
    savetf(TF2, fname1);
end

%% Third TF
if 1
    %% Third TF
    desc3 = sprintf('No bandpass');
    opts3 = tflab_options(1);
        opts3.tflab.loglevel = 1;
        opts3.td.window.width = pps;
        opts3.td.window.shift = pps;
        opts3.filestr = sprintf('%s-tf3',filestr);

    TF3 = tflab(Bo(:,1:2),Eo,opts3);
    % Modify default description of run
    TF3.Options.description = desc3;
    TF3.Metadata = meta;

    fname3 = fullfile(scriptdir(),'data','Middelpos',[opts3.filestr,'.mat']);
    savetf(TF3, fname3);
end

Middelpos_plot;
