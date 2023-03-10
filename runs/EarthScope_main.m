addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';

% Get input/output data
[B,E,t,infile,outfile] = EarthScope_clean(id);

%% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

%% Set output file base name using start/stop times of input data
filestr = sprintf('%s-%s-%s',id,...
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
iopts.info.chainid   = 'EarthScope';
iopts.info.stationid = id;


%% Compute first TF
desc1 = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts1 = tflab_options(1,iopts);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;
    opts1.filestr = sprintf('%s-tf1',filestr);

TF1 = tflab(B(:,1:2),E,opts1);
% Modify default description of run
TF1.Options.description = desc1;
TF1 = tflab_uncertainty(TF1);

fname1 = fullfile(scriptdir(),'data','EarthScope',id,[opts1.filestr,'.mat']);
savetf(TF1, fname1);
