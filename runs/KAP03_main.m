clear;

% Add location of tflab_setpaths() to path.
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
% Set all tflab paths.
tflab_setpaths(); 

chainid   = 'KAP03';
stationid = 'KAP103';

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP103_clean();

% Base file name for saved results
filestr = sprintf('%s-%s-%s',...
            stationid,datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));

% Set variable information (used for plots)
iopts = struct('info',struct(),'td',struct());
    iopts.info.instr     = {'$B_x$','$B_y$'};
    iopts.info.inunit    = 'nT';
    iopts.info.outstr    = {'$E_x$','$E_y$'};
    iopts.info.outunit   = 'mV/km';
    iopts.info.timeunit  = 's';
    iopts.info.timedelta = timedelta;
    iopts.info.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    iopts.info.chainid   = chainid;
    iopts.info.stationid = stationid;

logmsg('Interpolating over NaNs in E\n');
E = naninterp1(E);
logmsg('Interpolating over NaNs in B\n');
B = naninterp1(B);

% Make length an integer number of segments.
pps = 86400; % Points per segment
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

%% Compute transfer functions

%% First TF
desc1 = sprintf('OLS; %d %d-day segments',size(B,1)/pps,pps/86400);
opts1 = tflab_options(1,iopts);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = pps;
    opts1.td.window.shift = pps;
    opts1.filestr = sprintf('%s-tf1',filestr);

% Execute run
TF1 = tflab(B(:,1:2),E,opts1);

% Modify default description of run
TF1.Options.description = desc1;

% Compute uncertainties
TF1 = tflab_uncertainty(TF1);

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

    TF2 = read_edi([zread_dir,'/data/kap003.edi']);
    TF2.Z(TF2.Z > 1e31) = NaN;
    
    % Set options and data needed for metrics and plotting
    TF2.Options.info = TF1.Options.info;
    TF2.Options.filestr = sprintf('%s-tf2',filestr);
    TF2.Options.description = 'BIRP';
    TF2.In  = TF1.In;
    TF2.Out = TF1.Out;
    TF2.Time = TF1.Time;
    TF2 = tflab_metrics(TF2,TF1.Options,TF1.Segment.IndexRange,1);
end

fname = fullfile(scriptdir(),'data',chainid,stationid,[TF2.Options.filestr,'.mat']);
savetf(TF2,fname);

KAP03_plot;
