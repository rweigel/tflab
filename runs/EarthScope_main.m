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
    %edifile = 'ORF03bc_G3x.xml';
    Ikeep = [86400*12+1:86400*16];
    edifile = 'USArray.ORF03.2007.xml'; % Older version
    % http://ds.iris.edu/spud/emtf/14866915
    % http://ds.iris.edu/spudservice/data/14866913
end

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

if ~isempty(Ikeep)
    B = B(Ikeep,:);
    E = E(Ikeep,:);
    t = t(Ikeep);
end
%[B,E] = removemean(B,E);

E = bandpass_(E,[1/86400,0.5]);
B = bandpass_(B,[1/86400,0.5]);

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

%% First TF
tfn = 1;
pps = size(B,1);
desc = sprintf('OLS; 1 %d-day segment',pps/86400);
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.td.window.width = NaN;
    opts{tfn}.td.window.shift = NaN;
    %opts{tfn}.td.detrend.function = @bandpass_;
    %opts{tfn}.td.detrend.functionstr = 'bandpass_';
    opts{tfn}.td.detrend.functionargs = {[1/86400,0.5]};
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts{tfn}.description = desc;
    if short_run
        opts{tfn}.fd.bootstrap.N = NaN;
    end

TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});

TFs{tfn}.Metadata = meta; % Attach metadata used in plots

savetf(TFs{tfn}, fullfile(outdir, opts{tfn}.filestr));

%% Second TF
tfn = tfn + 1;
pps = 86400;
desc = sprintf('OLS; %d %d-day segment%s',size(B,1)/pps,pps/86400,plural(Ns));
opts{tfn} = tflab_options(1);
    opts{tfn}.tflab.loglevel = 1;
    opts{tfn}.td.window.width = pps;
    opts{tfn}.td.window.shift = pps;
    opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    opts{tfn}.description = desc;
    if short_run
        opts{tfn}.fd.bootstrap.N = NaN;
    end
TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});

TFs{tfn}.Metadata = meta; % Attach metadata used in plots
savetf(TFs{tfn}, fullfile(outdir, opts{tfn}.filestr));

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

    edifilefull = fullfile(scriptdir(),'data','EarthScope',id,'edi',edifile);
    EDI = read_edixml(edifilefull);

    % Set options and data needed for metrics and plotting
    TFs{tfn}.Metadata = meta;
    TFs{tfn}.Metadata.EDI = EDI;
    TFs{tfn}.Metadata.timestart = TFs{1}.Metadata.timestart;
    %TF3.Metadata.timestart = [TF3.Metadata.EDI.Start,'.000'];
    TFs{tfn}.Options.filestr = sprintf('%s-tf%d',filestr,tfn);

    % EMTF is listed as software at http://ds.iris.edu/spud/emtf/15014571
    TFs{tfn}.Options.description = 'EMTF';
    TFs{tfn}.Options.fd = opts{1}.fd;
    TFs{tfn}.Options.tflab.loglevel = opts{1}.tflab.loglevel;
    TFs{tfn}.Z = -transpose(TFs{tfn}.Metadata.EDI.Z); % Negative due to e^{+iwt} convention
    TFs{tfn}.ZVAR = transpose(TFs{tfn}.Metadata.EDI.ZVAR);
    TFs{tfn}.fe = 1./transpose(TFs{tfn}.Metadata.EDI.PERIOD);
    TFs{tfn}.In  = TFs{1}.In(:,1:2);
    TFs{tfn}.Out = TFs{1}.Out;

    TFs{tfn} = tflab_preprocess(TFs{tfn});
    TFs{tfn} = tflab_metrics(TFs{tfn});

    fname = fullfile(outdir,[TFs{tfn}.Options.filestr,'.mat']);
    savetf(TFs{tfn},fname);
end