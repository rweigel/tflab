addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

short_run = 1;

if 1
    id = 'VAQ58';
    edifile = 'VAQ58bc_FRDcoh.xml';
    Ikeep = [];
    %Ikeep = 1:4*86400;
    Ikeep = 1:6*86400;
    % http://ds.iris.edu/spud/emtf/15014571
    % http://ds.iris.edu/spudservice/data/15014570
end

if 1
    id = 'ORF03';
    %edifile = 'ORF03bc_G3x.xml'; % Older version
    %start = '2007-08-19T01:48:36';
    %stop = '2007-09-07T17:18:40';
    start = '2007-08-31T01:48:36';
    stop = '2007-09-04T01:48:35';
    edifile = 'USArray.ORF03.2007.xml'; % Name in XML linked to below.
    ediurl  = 'http://ds.iris.edu/spudservice/data/21636090';
    % http://ds.iris.edu/spud/emtf/14866915
end

if 0
    id = 'ORG03';
    start = '2007-08-31T01:48:36';
    stop = '2007-09-04T01:48:35';
    edifile = 'USArray.ORG03.2007.xml'; % Name in XML linked to below.
    ediurl = 'http://ds.iris.edu/spudservice/data/21636974';
    % http://ds.iris.edu/spud/emtf/14867137
end

dno = datenum(start,'yyyy-mm-ddTHH:MM:SS');
dnf = datenum(stop,'yyyy-mm-ddTHH:MM:SS');
dirstr  = sprintf('tfs-%s-%s',...
            datestr(dno,'yyyymmddTHHMMSS'),datestr(dnf,'yyyymmddTHHMMSS'));
rundir = fullfile(scriptdir(),'data','EarthScope',id,dirstr);

% Get input/output data
[B,E,t,infile,outfile] = EarthScope_clean(id);

Ik = t >= dno & t <= dnf;
B = B(Ik,:);
E = E(Ik,:);
t = t(Ik);

[B,E] = removemean(B,E);

E = bandpass_(E,[1/86400,0.5]);
B = bandpass_(B,[1/86400,0.5]);

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
    meta.timestop  = datestr(t(end),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'EarthScope';
    meta.stationid = id;

const_term = 0;

%% First TF
tfn = 1;
pps = size(B,1);
desc = sprintf('OLS; 1 %d-day segment',pps/86400);
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.td.window.width = NaN;
    opts{tfn}.td.window.shift = NaN;
    opts{tfn}.fd.regression.const_term = const_term;
    %opts{tfn}.td.detrend.function = @bandpass_;
    %opts{tfn}.td.detrend.functionstr = 'bandpass_';
    opts{tfn}.td.detrend.functionargs = {[1/86400,0.5]};
    opts{tfn}.filestr = sprintf('%s-tf%d',id,tfn);
    opts{tfn}.description = desc;
    if short_run
        opts{tfn}.fd.bootstrap.N = NaN;
    end

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});

TFs{tfn}.Metadata = meta; % Attach metadata used in plots

savetf(TFs{tfn}, fullfile(rundir, opts{tfn}.filestr));

%% Second TF
tfn = tfn + 1;
pps = 86400;
desc = sprintf('OLS; %d %d-day segment%s',size(B,1)/pps,pps/86400,plural(pps/86400));
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.fd.regression.const_term = const_term;
    opts{tfn}.filestr = sprintf('%s-tf%d',id,tfn);
    opts{tfn}.description = desc;
    if short_run
        opts{tfn}.fd.bootstrap.N = NaN;
    end
TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
TFs{tfn}.Metadata = meta; % Attach metadata used in plots
savetf(TFs{tfn}, fullfile(rundir, opts{tfn}.filestr));

if 1
    tfn = tfn + 1;
    %% TF computed using EMTF
    zread_dir = fullfile(scriptdir(),'zread');
    if ~exist(zread_dir,'dir')
        url = 'https://github.com/rweigel/zread';
        cmd = 'cd %s; git clone %s; cd zread; git checkout master; git checkout 1519ed4';
        cmd = sprintf(cmd,scriptdir(),url);
        fprintf('Calling system with command: %s\n',cmd);
        [status,msg] = system(cmd);
        if status ~= 0
            fprintf(['Command failed. Download and unzip '...
                    'github.com/rweigel/zread/archive/refs/heads/master.zip'...
                    'in this directory\n']);
            error('System command failed: %s\nMessage:\n%s\n',cmd,msg);
        end
    end
    addpath(zread_dir);

    edifile = fullfile('data','EarthScope',id,'edi',edifile);
    edifilefull = fullfile(scriptdir(),edifile);
    if ~exist(edifilefull,'file')
        if ~exist(fileparts(edifilefull),'dir')
            mkdir(fileparts(edifilefull));
        end
        fprintf('Downloading %s to %s\n',ediurl,edifilefull);
        websave(edifilefull, ediurl);
    end
    EDI = read_edixml(edifilefull);

    % Set options and data needed for metrics and plotting
    TFs{tfn}.Metadata = meta;
    TFs{tfn}.Metadata.EDI = EDI;
    TFs{tfn}.Metadata.timestart = TFs{1}.Metadata.timestart;
    %TF3.Metadata.timestart = [TF3.Metadata.EDI.Start,'.000'];
    TFs{tfn}.Options.filestr = sprintf('%s-tf%d',id,tfn);

    % EMTF is listed as software at http://ds.iris.edu/spud/emtf/15014571
    TFs{tfn}.Options.description = 'EMTF';
    TFs{tfn}.Options.fd = opts{1}.fd;
    TFs{tfn}.Options.tflab.loglevel = opts{1}.tflab.loglevel;
    TFs{tfn}.Z = -transpose(TFs{tfn}.Metadata.EDI.Z); % Negative due to e^{+iwt} convention
    TFs{tfn}.ZVAR = transpose(TFs{tfn}.Metadata.EDI.ZVAR);
    TFs{tfn}.fe = 1./transpose(TFs{tfn}.Metadata.EDI.PERIOD);
    TFs{tfn}.In  = TFs{1}.In(:,1:2);
    TFs{tfn}.Out = TFs{1}.Out;

    if const_term
        tmp = nan(size(TFs{tfn}.Z,1),1);
        TFs{tfn}.Z = [TFs{tfn}.Z(:,1:2),tmp,TFs{tfn}.Z(:,3:4),tmp];
    end

    TFs{tfn} = tflab_preprocess(TFs{tfn});
    TFs{tfn} = tflab_metrics(TFs{tfn});

    fname = fullfile(rundir,[TFs{tfn}.Options.filestr,'.mat']);
    savetf(TFs{tfn},fname);
end