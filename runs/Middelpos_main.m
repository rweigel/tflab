clear;
close all % To reduce memory.
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

short_run = 0;

if short_run
    Nboot = NaN;
    start = '20120712';
    stop = '20120716';
else
    Nboot = 100;
    start = '20120712';
    stop = '20121107';
end

%% Set output file base name using start/stop times of input data
filestr = 'Middelpos';
dirstr  = sprintf('tfs-%s-%s',start,stop);
rundir = fullfile(scriptdir(),'data',filestr,dirstr);

% Read input/output data
[B,E,t] = Middelpos_clean(start,stop,rundir);

% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

if short_run
    B = B(1:pps*5,:);
    E = E(1:pps*5,:);
    t = t(1:pps*5);
end

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
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

%% TF1
tfn = 1;
logmsg('-- Computing TF%d --\n',tfn);
pps = 86400;
desc = sprintf('OLS; One %d-day segment',size(B,1)/pps);
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
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
desc = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts{tfn} = tflab_options(1);
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.description = desc;
    opts{tfn}.fd.bootstrap.N = Nboot;
    opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
TFs{tfn} = [];

%% TF3
tfn = 3;
logmsg('-- Computing TF%d --\n',tfn);
Tm = 1*86400;
band = [1/Tm,0.5];
desc = 'OLS; One %d-day segment; %d-day bandpass;';
desc = sprintf(desc,size(B,1)/pps,Tm/86400);

opts{tfn} = tflab_options(1);
    opts{tfn}.td.detrend.function = @bandpass_;
    opts{tfn}.td.detrend.functionstr = 'bandpass_';
    opts{tfn}.td.detrend.functionargs = {band};
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
desc = 'OLS; %d %d-day segments; %d-day bandpass';
desc = sprintf(desc,size(B,1)/pps,pps/86400,Tm/86400);

opts{tfn} = tflab_options(1);
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.td.detrend.function = @bandpass_;
    opts{tfn}.td.detrend.functionstr = 'bandpass_';
    opts{tfn}.td.detrend.functionargs = {band};
    opts{tfn}.description = desc;
    opts{tfn}.fd.bootstrap.N = Nboot;
    opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
TFs{tfn} = [];