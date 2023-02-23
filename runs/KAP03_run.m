clear;

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

stationid = 'KAP103';
chainid = 'KAP03';

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP103_clean();

filestr = sprintf('%s-%s-%s',...
            stationid,datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));

% Variable name information
iopts = struct('info',struct(),'td',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/km';
iopts.info.timeunit = 's';
iopts.info.timedelta = timedelta;
iopts.info.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
iopts.info.stationid = stationid;
iopts.info.chainid = chainid;

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
S1 = tflab(B(:,1:2),E,opts1);

% Modify default description of run
S1.Options.description = desc1;

% Compute uncertainties
S1 = tflab_uncertainty(S1);

% Save results
fname = fullfile(scriptdir(),'data',stationid,[S1.Options.filestr,'.mat']);
savetf(S1,fname);

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

    S2 = read_edi([zread_dir,'/data/kap003.edi']);
    S2.Z(S2.Z > 1e31) = NaN;
    
    % Set options and data needed for metrics and plotting
    S2.Options.info = S1.Options.info;
    S2.Options.filestr = sprintf('%s-tf2',filestr);
    S2.Options.description = 'BIRP';
    S2.In  = S1.In;
    S2.Out = S1.Out;
    S2.Time = S1.Time;
    S2 = tflab_metrics(S2,S1.Options,S1.Segment.IndexRange,1);
end

fname = fullfile(scriptdir(),'data',stationid,[S2.Options.filestr,'.mat']);
savetf(S2,fname);
