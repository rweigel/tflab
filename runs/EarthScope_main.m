addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';
edifile = 'VAQ58bc_FRDcoh.xml';
% http://ds.iris.edu/spud/emtf/15014571
% http://ds.iris.edu/spudservice/data/15014570

%id = 'ORF03';
%edifile = 'ORF03bc_G3x.xml';
%edifile = 'USArray.ORF03.2007.xml'; % Older version
% http://ds.iris.edu/spud/emtf/14866915
% http://ds.iris.edu/spudservice/data/14866913

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

if strcmp(id,'VAQ58')
    B = B(1:6*86400,:);
    E = E(1:6*86400,:);
    t = t(1:6*86400);
end

if 1 || strcmp(id,'VAQ58')
    %% Band pass
    Tm = 2*86400;
    band = [1/Tm,0.5];
    B = bandpass_(B,band);
    E = bandpass_(E,band);
    %keyboard
    %E = E(Tm+1:end-Tm,:);
    %B = B(Tm+1:end-Tm,:);
    %t = t(Tm+1:end-Tm,:);
end

if 0
    I = find(t > 736498.3 & t < 736502);
    B = B(I,:);
    E = E(I,:);
    t = t(I);
end

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
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'EarthScope';
    meta.stationid = id;

if 1
%% Compute first TF
tfn = 1;
Ns = size(B,1)/pps;
%desc1 = sprintf('TFLab; %d %d-day segment%s',Ns,pps/86400,plural(Ns));
desc1 = sprintf('TFLab');
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;
    opts1.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts1.description = desc1;
TF1 = tflab(B(:,1:2),E,opts1);

TF1.Metadata = meta; % Attach metadata used in plots

savetf(TF1, fullfile(outdir, opts1.filestr));
end


%% Compute second TF
tfn = 2;
Ns = 1;
pps = size(B,1);
%desc2 = sprintf('TFLab; %d %d-day segment%s',Ns,pps/86400,plural(Ns));
desc2 = sprintf('TFLab');
opts2 = tflab_options(1);
    opts2.tflab.loglevel = 1;
    opts2.td.window.width = NaN;
    opts2.td.window.shift = NaN;
    opts2.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts2.description = desc2;
TF2 = tflab(B(:,1:2),E,opts2);

TF2.Metadata = meta; % Attach metadata used in plots

savetf(TF2, fullfile(outdir, opts2.filestr));

%zplot(TF2)

if 1
%% Read TF computed using EMTF
zread_dir = [fileparts(mfilename('fullpath')),'/zread'];
if ~exist(zread_dir,'dir')
    url = 'https://github.com/rweigel/zread';
    com = sprintf('cd %s; git clone %s; cd zread; git checkout master; git checkout 1519ed4',...
                  fileparts(mfilename('fullpath')),url);
    fprintf('Calling system with command: %s\n',com);
    [status,msg] = system(com);
    if status ~= 0
        fprintf(['Command failed. Download and unzip '...
                 'github.com/rweigel/zread/archive/refs/heads/master.zip'...
                 'in this directory\n']);
        error('System command failed: %s\nMessage:\n%s\n',com,msg);            
    end
end
addpath(zread_dir);

edifilefull = fullfile(scriptdir(),'data','EarthScope',id,'edi',edifile);
EDI = read_edixml(edifilefull);

% Set options and data needed for metrics and plotting
TF3.Metadata = meta;
TF3.Metadata.EDI = EDI;
TF3.Metadata.timestart = TF2.Metadata.timestart;
%TF3.Metadata.timestart = [TF3.Metadata.EDI.Start,'.000'];
TF3.Options.filestr = sprintf('%s-tf3',filestr);

% EMTF is listed as software at http://ds.iris.edu/spud/emtf/15014571
TF3.Options.description = 'EMTF';
TF3.Options.fd = opts1.fd;
TF3.Options.tflab.loglevel = opts1.tflab.loglevel;
TF3.Z = -transpose(TF3.Metadata.EDI.Z); % Negative due to e^{+iwt} convention
TF3.ZVAR = transpose(TF3.Metadata.EDI.ZVAR);
TF3.fe = 1./transpose(TF3.Metadata.EDI.PERIOD);
TF3.In  = TF1.In(:,1:2); 
TF3.Out = TF1.Out;

TF3 = tflab_metrics(TF3);

fname = fullfile(outdir,[TF3.Options.filestr,'.mat']);
savetf(TF3,fname);
end