function Middelpos_main(rundir, filestr, start, stop, const_term, Nboot)

% Read input/output data
[B,E,t] = Middelpos_clean(start,stop,rundir,0);

% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

%% Set common metadata
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = 1;
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.timestop  = datestr(t(end),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

%% TF1
tfn = 1;
logmsg('-- Computing TF%d --\n',tfn);
pps = 86400;
desc = sprintf('1 %d-day segment',size(B,1)/pps);
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.fd.regression.const_term = const_term;
    opts{tfn}.description = desc;
    opts{tfn}.fd.bootstrap.N = Nboot;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
TFs{tfn} = [];

%% TF2
tfn = 2;
logmsg('-- Computing TF%d --\n',tfn);
pps = 86400;
desc = sprintf('%d %d-day segments',size(B,1)/pps,pps/86400);
opts{tfn} = tflab_options(1);
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.fd.regression.const_term = const_term;
    opts{tfn}.description = desc;
    opts{tfn}.fd.bootstrap.N = Nboot;
    opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
TFs{tfn} = [];

if 1
    %% TF3
    tfn = 3;
    logmsg('-- Computing TF%d --\n',tfn);
    Tm = 1*86400;
    band = [1/Tm,0.5];
    desc = '1 %d-day seg.; %d-day bandpass;';
    desc = sprintf(desc,size(B,1)/pps,Tm/86400);
    
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.detrend.function = @bandpass_;
        opts{tfn}.td.detrend.functionstr = 'bandpass_';
        opts{tfn}.td.detrend.functionargs = {band};
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
    
    %% TF4
    tfn = 4;
    logmsg('-- Computing TF%d --\n',tfn);
    pps = 86400;
    Tm = 1*86400;
    band = [1/Tm,0.5];
    desc = '%d %d-day segments; %d-day bandpass';
    desc = sprintf(desc,size(B,1)/pps,pps/86400,Tm/86400);
    
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.window.width = pps;
        opts{tfn}.td.window.shift = pps;
        opts{tfn}.td.detrend.function = @bandpass_;
        opts{tfn}.td.detrend.functionstr = 'bandpass_';
        opts{tfn}.td.detrend.functionargs = {band};
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
end

if 0    
    %% TF5
    tfn = 5;
    logmsg('-- Computing TF%d --\n',tfn);
    meta.instr     = {'$B_x$'};
    meta.outstr    = {'$E_y$'};
    pps = 86400;
    desc = sprintf('%d %d-day segments',size(B,1)/pps,pps/86400);
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.window.width = pps;
        opts{tfn}.td.window.shift = pps;
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    TFs{tfn} = tflab(B(:,1),E(:,2),opts{tfn});
    TFs{tfn}.Metadata = meta;
    
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
end
