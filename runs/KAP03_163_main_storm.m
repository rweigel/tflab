clear;

chainid   = 'KAP03';
stationid = 'KAP163';

KAP03_main;

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP03_clean(stationid);

Nr = (86400/5)*3; % Number of records
I1a = [1:Nr]; % First 3 days
B1a = B(I1a,:);
E1a = E(I1a,:);
t1a = t(I1a);

filestr1a = sprintf('%s-%s-%s',...
            stationid,datestr(t1a(1),'yyyymmdd'),...
            datestr(t1a(end),'yyyymmdd'));

I1b = [Nr+1:2*Nr]; % Second Nd days
B1b = B(I1b,:);
E1b = E(I1b,:);
t1b = t(I1b);

% Base file name for saved results
filestr1b = sprintf('%s-%s-%s',...
            stationid,datestr(t1b(1),'yyyymmdd'),...
            datestr(t1b(end),'yyyymmdd'));

%% Compute transfer functions

%% First TF
meta.timestart = datestr(t1a(1),'yyyy-mm-ddTHH:MM:SS.FFF');
opts1a = tflab_options(1);
    opts1a.tflab.loglevel = 1;
    opts1a.td.window.width = NaN;
    opts1a.td.window.shift = NaN;
    opts1a.description = 'TFLab/3-day storm';
    opts1a.filestr = sprintf('%s-tf1a',filestr1a);
    
% Execute run
TF1a = tflab(B1a(:,1:2),E1a,opts1a);
TF1a.Metadata = meta;

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF1a.Options.filestr,'.mat']);
savetf(TF1a,fname);

%% Second TF
meta.timestart = datestr(t1b(1),'yyyy-mm-ddTHH:MM:SS.FFF');
opts1b = tflab_options(1);
    opts1b.tflab.loglevel = 1;
    opts1b.td.window.width = NaN;
    opts1b.td.window.shift = NaN;
    opts1b.description = 'TFLab/3-day quiet';
    opts1b.filestr = sprintf('%s-tf1b',filestr1b);
    
% Execute run
TF1b = tflab(B1b(:,1:2),E1b,opts1b);
TF1b.Metadata = meta;

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

    EDI = read_edi([scriptdir(),'/data/KAP03/edi/',lower(stationid),'.edi']);
    
    TF2.Metadata = meta;
    TF2.Metadata.EDI = EDI;
    TF2.Z = TF2.Metadata.EDI.Z;
    TF2.Z(TF2.Z > 1e31) = NaN;
    TF2.ZVAR = TF2.Metadata.EDI.ZVAR;

    TF2a = TF2;
    TF2b = TF2;
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
    TF2a.Options.filestr = sprintf('%s-tf2a',filestr);
    TF2b.Options.filestr = sprintf('%s-tf2b',filestr);
    TF2a.Options.description = 'BIRP/All/3-day storm';
    TF2b.Options.description = 'BIRP/All/3-day quiet';
    TF2a.Options.fd = opts1a.fd;
    TF2b.Options.fd = opts1a.fd;
    TF2a.Options.tflab.loglevel = opts1a.tflab.loglevel;
    TF2b.Options.tflab.loglevel = opts1a.tflab.loglevel;
    % Convert frequencies to be consistent with cadence of time series
    % New frequency has units of 1/sample
    TF2a.fe = TF2.Metadata.EDI.fe*timedelta;
    TF2b.fe = TF2.Metadata.EDI.fe*timedelta;
    TF2a.In  = B(I1a(1):I1a(end),1:2);
    TF2b.In  = B(I1b(1):I1b(end),1:2);
    TF2a.Out = E(I1a(1):I1a(end),:);
    TF2b.Out = E(I1b(1):I1b(end),:);
    %TF2.Z = TF2.Z;
    TF2a = tflab_metrics(TF2a);
    TF2b = tflab_metrics(TF2b);
end

fnamea = fullfile(scriptdir(),'data',chainid,stationid,[TF2a.Options.filestr,'.mat']);
savetf(TF2a,fnamea);
fnameb = fullfile(scriptdir(),'data',chainid,stationid,[TF2b.Options.filestr,'.mat']);
savetf(TF2b,fnameb);

KAP03_plot;
