addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';
% http://ds.iris.edu/spud/emtf/15014571
% http://ds.iris.edu/spudservice/data/15014570

outdir = fullfile(scriptdir(),'data','EarthScope',id);

% Get input/output data
[B,E,t,infile,outfile] = EarthScope_clean(id);

%% Make length an integer number of segments.
pps = 86400;                % Points per segment
Ns  = floor(size(B,1)/pps); % Number of segments
I = pps*Ns;
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

B = B(1:6*86400,:);
E = E(1:6*86400,:);
t = t(1:6*86400);

%% Set output file base name using start/stop times of input data
filestr = sprintf('%s-%s-%s',id,...
                    datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));

%% Set metadata used for plots
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = 1;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'EarthScope';
    meta.stationid = id;

%% Compute first TF
tfn = 1;
Ns = size(B,1)/pps;
desc1 = sprintf('OLS; %d %d-day segment%s',Ns,pps/86400,plural(Ns));
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;
    opts1.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts1.description = desc1;
TF1 = tflab(B(:,1:2),E,opts1);

%TF1 = tflab_uncertainty(TF1);

TF1.Metadata = meta; % Attach metadata used in plots

savetf(TF1, fullfile(outdir, opts1.filestr));

%% Compute second TF
Ns = 1;
pps = size(B,1);
desc2 = sprintf('OLS; %d %d-day segment%s',Ns,pps/86400,plural(Ns));
opts2 = tflab_options(1);
    opts2.tflab.loglevel = 1;
    opts2.td.window.width = pps;
    opts2.td.window.shift = pps;
    opts2.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts2.description = desc2;
TF2 = tflab(B(:,1:2),E,opts2);

TF2.Metadata = meta; % Attach metadata used in plots

savetf(TF2, fullfile(outdir, opts2.filestr));
