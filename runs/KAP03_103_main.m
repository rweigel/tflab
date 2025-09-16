clear;

chainid   = 'KAP03';
stationid = 'KAP103';

KAP03_main;

ppd = 86400/5; % Points per day
if 0
    nd = floor(size(E,1)/ppd);
else
    nd = 10;
end
E = E(1:nd*ppd,:);
B = B(1:nd*ppd,:);
t = t(1:nd*ppd);

to = datestr(t(1),'yyyymmdd');
tf = datestr(t(end),'yyyymmdd');

outdir = fullfile(scriptdir(),'data',chainid,stationid,sprintf('tfs-%s-%s',to,tf));

%% First TF
opts1 = tflab_options(1);
    opts1.tflab.loglevel = 1;
    opts1.td.window.width = NaN;
    opts1.td.window.shift = NaN;
    opts1.description = sprintf('1 %d-day segment',nd);
    opts1.filestr = sprintf('%s-tf1',stationid);

% Execute run
TF1 = tflab(B(:,1:2),E,opts1);
TF1.Metadata = meta;

% Save results
fname = fullfile(outdir,[TF1.Options.filestr,'.mat']);
savetf(TF1,fname);

%% Second TF
opts2 = tflab_options(1);
    opts2.tflab.loglevel = 1;
    opts2.td.window.width = ppd;
    opts2.td.window.shift = ppd;
    opts2.description = sprintf('%d 1-day segment',nd);
    opts2.filestr = sprintf('%s-tf2',stationid);

% Execute run
TF2 = tflab(B(:,1:2),E,opts2);
TF2.Metadata = meta;

% Save results
fname = fullfile(outdir,[TF2.Options.filestr,'.mat']);
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
    TF3.Z(TF3.Z > 1e31) = NaN; % TODO: Get from EMPTY= line in file; kap003 file has EMPTY= 0.1000000E+33
    TF3.ZVAR = TF3.Metadata.EDI.ZVAR;
    TF3.ZVAR(TF3.ZVAR > 1e31) = NaN;

    % edi file read has filenam=kap163/kap163as.{ex,ey,...}
    % and the file we read is kap163as.ts. edi file also nas
    % nread = 463807 and kap163as.ts has 464969 lines, so it appears
    % edi file was based on about full range of data. 
    
    % Set options and data needed for metrics and plotting
    TF3.Options.filestr = sprintf('%s-tf3',stationid);
    TF3.Options.description = 'BIRRP';
    TF3.Options.fd.evalfreq = opts1.fd.evalfreq;
    TF3.Options.fd.window = opts1.fd.window;
    TF3.Options.fd.regression.const_term = 0;
    TF3.Options.tflab.loglevel = opts1.tflab.loglevel;
    % Convert frequencies to be consistent with cadence of time series
    % 
    TF3.fe = TF3.Metadata.EDI.fe*timedelta;
    TF3.In  = B(:,1:2);
    TF3.Out = E;
    %TF3 = tflab_metrics(TF3);
end

fnamea = fullfile(outdir,[TF3.Options.filestr,'.mat']);
savetf(TF3,fnamea);

%KAP03_plot;
