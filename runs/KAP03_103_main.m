clear;

chainid   = 'KAP03';
stationid = 'KAP103';

KAP03_main;

%% Compute transfer functions

%% First TF
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = NaN;
    opts1.td.window.shift = NaN;
    opts1.description = 'TFLab 1';
    opts1.filestr = sprintf('%s-%s-%s-tf1',...
                        stationid,datestr(t(1),'yyyymmdd'),...
                        datestr(t(end),'yyyymmdd'));    
% Execute run
TF1 = tflab(B(:,1:2),E,opts1);
TF1.Metadata = meta;

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF1.Options.filestr,'.mat']);
savetf(TF1,fname);

%% Second TF
opts2 = tflab_options(1);
    opts2.tflab.loglevel = 1;
    opts2.td.window.width = 86400/5;
    opts2.td.window.shift = 86400/5;
    opts2.description = 'TFLab 2';
    opts2.filestr = sprintf('%s-%s-%s-tf2',...
                        stationid,datestr(t(1),'yyyymmdd'),...
                        datestr(t(end),'yyyymmdd'));    
% Execute run
TF2 = tflab(B(:,1:2),E,opts2);
TF2.Metadata = meta;

% Save results
fname = fullfile(scriptdir(),'data',chainid,stationid,[TF2.Options.filestr,'.mat']);
savetf(TF2,fname);

%% Read TF computed using BIRRP
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
    
    TF3.Metadata = meta;
    TF3.Metadata.EDI = EDI;
    TF3.Z = TF3.Metadata.EDI.Z;
    TF3.Z(TF3.Z > 1e31) = NaN;
    TF3.ZVAR = TF3.Metadata.EDI.ZVAR;

    filestr = sprintf('%s-unknown-unknown',stationid);
    if strcmp(stationid,'KAP163')
        % edi file read has filenam=kap163/kap163as.{ex,ey,...}
        % and the file we read is kap163as.ts. edi file also nas
        % nread = 463807 and kap163as.ts has 464969 lines, so it appears
        % edi file was based on about full range of data. 
        filestr = sprintf('%s-%s-%s',...
                  stationid,datestr(t(1),'yyyymmdd'),...
                  datestr(t(end),'yyyymmdd'));
    end
    
    % Set options and data needed for metrics and plotting
    TF3.Options.filestr = sprintf('%s-tf3',filestr);
    TF3.Options.description = 'BIRRP';
    TF3.Options.fd = opts1.fd;
    TF3.Options.tflab.loglevel = opts1.tflab.loglevel;
    % Convert frequencies to be consistent with cadence of time series
    % New frequency has units of 1/sample
    TF3.fe = TF3.Metadata.EDI.fe*timedelta;
    TF3.In  = B(:,1:2);
    TF3.Out = E;
    TF3 = tflab_metrics(TF3);
end

fnamea = fullfile(scriptdir(),'data',chainid,stationid,[TF3.Options.filestr,'.mat']);
savetf(TF3,fnamea);

%KAP03_plot;
