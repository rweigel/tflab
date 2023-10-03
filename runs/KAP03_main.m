clear;

% Add location of tflab_setpaths() to path.
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
% Set all tflab paths.
tflab_setpaths(); 

chainid   = 'KAP03';
stationid = 'KAP163';

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP03_clean(stationid);

I1a = [1:(86400/5)]*3; % First 3 days
B1a = B(I1a,:);
E1a = E(I1a,:);
t1a = t(I1a);

filestr1a = sprintf('%s-%s-%s',...
            stationid,datestr(t1a(1),'yyyymmdd'),datestr(t1a(end),'yyyymmdd'));

I1b = [I1a(end)+1:I1a(end)+(86400/5)*3]; % Second 3 days
B1b = B(I1b,:);
E1b = E(I1b,:);
t1b = t(I1b);

% Base file name for saved results
filestr1b = sprintf('%s-%s-%s',...
            stationid,datestr(t1b(1),'yyyymmdd'),datestr(t1b(end),'yyyymmdd'));

% Set variable information (used for plots)
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = timedelta;
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1/timedelta; % Multiply frequencies by this to get frequnit
    meta.timestart = datestr(t1a(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = chainid;
    meta.stationid = stationid;


%% Compute transfer functions

%% First TF
opts1a = tflab_options(1);
    opts1a.tflab.loglevel = 1;
    opts1a.td.window.width = NaN;
    opts1a.td.window.shift = NaN;
    opts1a.filestr = sprintf('%s-tf1a',filestr1a);
    
% Execute run
TF1a = tflab(B1a(:,1:2),E1a,opts1a);
TF1a.Metadata = meta;

% Compute uncertainties
%TF1 = tflab_uncertainty(TF1);

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF1a.Options.filestr,'.mat']);
savetf(TF1a,fname);

%% Second TF
meta.timestart = datestr(t1b(1),'yyyy-mm-ddTHH:MM:SS.FFF');
opts1b = tflab_options(1);
    opts1b.tflab.loglevel = 1;
    opts1b.td.window.width = NaN;
    opts1b.td.window.shift = NaN;
    opts1b.filestr = sprintf('%s-tf1b',filestr1b);
    
% Execute run
TF1b = tflab(B1b(:,1:2),E1b,opts1b);
TF1b.Metadata = meta;

% Compute uncertainties
%TF1 = tflab_uncertainty(TF1);

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF1b.Options.filestr,'.mat']);
savetf(TF1b,fname);

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

    filestr = sprintf('%s-unknown-unknown',stationid);
    if strcmp(stationid,'KAP163')
        % edi file read has filenam=kap163/kap163as.{ex,ey,...}
        % and the file we read is kap163as.ts. edi file also nas
        % nread = 463807 and kap163as.ts has 464969 lines, so it appears
        % edi file was based on about full range of data. 
        filestr = sprintf('%s-%s-%s',...
                  stationid,datestr(t(1),'yyyymmdd'),datestr(t(end),'yyyymmdd'));
    end
    
    % Set options and data needed for metrics and plotting
    TF2.Metadata = meta;
    TF2.Metadata.freqsf = 1;
    TF2.Options.filestr = sprintf('%s-tf2',filestr);
    TF2.Options.description = 'BIRP';
    TF2.Options.fd = opts1a.fd;
    TF2.Options.tflab.loglevel = opts1a.tflab.loglevel;
    TF2.In  = TF1a.In;
    TF2.Out = TF1a.Out;
    TF2.Z = TF2.Z;
    TF2 = tflab_metrics(TF2);
end

fname = fullfile(scriptdir(),'data',chainid,stationid,[TF2.Options.filestr,'.mat']);
savetf(TF2,fname);

KAP03_plot;
