clear;

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

% Get input/output data
[B,E,t,infile,outfile] = Middelpos_clean();

%% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

if 1
    % Do short segment for testing
    B = B(1:pps*5,:);
    E = E(1:pps*5,:);
    t = t(1:pps*5);
end

%% Set output file base name using start/stop times of input data
fmt = 'yyyymmdd';
filestr = 'Middelpos';
dirstr  = sprintf('tfs-%s-%s',datestr(t(1),fmt),datestr(t(end),fmt));
basepath = fullfile(scriptdir(),'data',filestr,dirstr);

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
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

%% TF1
tfn = 1;
logmsg('Computing TF%d\n',tfn);
pps = 86400;
desc = sprintf('OLS; One %d-day segment',size(B,1)/pps);
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.description = desc;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(basepath,opts{tfn}.filestr));

%% TF2
tfn = 2;
logmsg('Computing TF%d\n',tfn);
pps = 86400;
desc = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts{tfn} = tflab_options(1);
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.description = desc;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(basepath,opts{tfn}.filestr));

%% TF3
tfn = 3;
logmsg('Computing TF%d\n',tfn);
Tm = 1*86400;
band = [1/Tm,0.5];
B_ = bandpass_(B,band);
E_ = bandpass_(E,band);
%E_ = E(Tm+1:end-Tm,:);
%B_ = B(Tm+1:end-Tm,:);
desc = 'OLS; One %d-day segment; %d-day bandpass;';
desc = sprintf(desc,size(B,1)/pps,Tm/86400);

opts{tfn} = tflab_options(1);
    opts{tfn}.description = desc;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B_(:,1:2),E_,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(basepath,opts{tfn}.filestr));

%% TF4
tfn = 4;
logmsg('Computing TF%d\n',tfn);
pps = 86400;
Tm = 1*86400;
band = [1/Tm,0.5];
B_ = bandpass_(B,band);
E_ = bandpass_(E,band);
%E_ = E(Tm+1:end-Tm,:);
%B_ = B(Tm+1:end-Tm,:);
desc = 'OLS; %d %d-day segments; %d-day bandpass';
desc = sprintf(desc,size(B,1)/pps,pps/86400,Tm/86400);

opts{tfn} = tflab_options(1);
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.description = desc;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);

TFs{tfn} = tflab(B_(:,1:2),E_,opts{tfn});
TFs{tfn}.Metadata = meta;

savetf(TFs{tfn}, fullfile(basepath,opts{tfn}.filestr));
