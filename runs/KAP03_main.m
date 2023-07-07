clear;

% Add location of tflab_setpaths() to path.
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
% Set all tflab paths.
tflab_setpaths(); 

chainid   = 'KAP03';
stationid = 'KAP163';

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP03_clean(stationid);

if 0
    %% Band pass
    addpath(fullfile(scriptdir(),'..','fft'));
    Tm = (86400/5)/4;
    band = [1/Tm,0.5];
    B = bandpass_(B,band);
    E = bandpass_(E,band);
    %Eb = E(Tm+1:end-Tm,:);
    %Bb = B(Tm+1:end-Tm,:);
end

% Make length an integer number of segments.
pps = 86400; % Points per segment
%I = pps*floor(size(B,1)/pps);
I = [1:86400];
B = B(I,:);
E = E(I,:);
t = t(I);

% Base file name for saved results
filestr = sprintf('%s-%s-%s',...
            stationid,datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));

% Set variable information (used for plots)
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = timedelta;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = chainid;
    meta.stationid = stationid;


%% Compute transfer functions

%% First TF
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = NaN;
    opts1.td.window.shift = NaN;
    opts1.filestr = sprintf('%s-tf1',filestr);
    
% Execute run
TF1 = tflab(B(:,1:2),E,opts1);
TF1.Metadata = meta;

% Compute uncertainties
%TF1 = tflab_uncertainty(TF1);

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF1.Options.filestr,'.mat']);
savetf(TF1,fname);

%% Read TF computed using BIRP
if strcmp(chainid,'KAP03')
    zread_dir = [fileparts(mfilename('fullpath')),'/zread'];
    if ~exist(zread_dir,'dir')
        url = 'https://github.com/rweigel/zread';
        com = sprintf('cd %s; git clone %s; cd zread; git checkout 8be764e50439db308bbb0b51b886bf0b7fb10c24',fileparts(mfilename('fullpath')),url);
        fprintf('Calling system with command: %s\n',com);
        [status,msg] = system(com);
        if status ~= 0
            fprintf('Command failed. Download and unzip https://github.com/rweigel/zread/archive/refs/heads/master.zip in this directory\n');
            error('System command failed: %s\nMessage:\n%s\n',com,msg);            
        end
    end
    addpath(zread_dir);

    TF2 = read_edi([scriptdir(),'/data/KAP03/edi/',lower(stationid),'.edi']);
    TF2.Z(TF2.Z > 1e31) = NaN;
    
    % Set options and data needed for metrics and plotting
    TF2.Metadata = meta;
    TF2.Options.filestr = sprintf('%s-tf2',filestr);
    TF2.Options.description = 'BIRP';
    TF2.Options.fd = opts1.fd;
    TF2.Options.tflab.loglevel = opts1.tflab.loglevel;
    TF2.In  = TF1.In;
    TF2.Out = TF1.Out;
    TF2 = tflab_metrics(TF2);
end

fname = fullfile(scriptdir(),'data',chainid,stationid,[TF2.Options.filestr,'.mat']);
savetf(TF2,fname);

KAP03_plot;
