clear;

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

chainid = 'KAP03';
print_figs = 0;

if 0
    siteid = 'KAP163';
    tfs = {...
            '20031028-20031124-tf1',...
            '20031028-20031124-tf2'...
            '20031028-20031124-tf3'...
           };
end

if 0
    siteid = 'KAP163';
    tfs = {...
            '20031028-20031031-tf1a',...
            '20031031-20031103-tf1b',...
            '20031028-20031124-tf2a',...
            '20031028-20031124-tf2b'...
           };
    tfs = tfs(1:3);
end

if 0
    siteid = 'KAP103';
    rundir = 'tfs-20031108-20031118';
    tfs = {'tf1','tf2','tf3'};
end

if 1
    siteid = 'KAP103';
    rundir = 'tfs-20031108-20031204';
    tfs = {'tf1','tf2','tf3'};
end

if 0
    siteid = 'KAP103';
    rundir = 'tfs-20031121-20031204';
    tfs = {'tf1','tf2','tf3'};
end

outdir = fullfile(scriptdir(),'data','KAP03',siteid,rundir);
for tfn = 1:length(tfs)
    f = fullfile(outdir, sprintf('%s-%s.mat',siteid,tfs{tfn}));
    TFs{tfn} = loadtf(f);
    TFs{tfn} = tflab_preprocess(TFs{tfn});
    TFs{tfn} = tflab_metrics(TFs{tfn});
end

%% Common plot options
copts.print = print_figs;  % Set to 1 to print pdf of each figure created.
copts.printOptions.printDir = fullfile(outdir,'figures');
copts.printOptions.printFormats = {'pdf','png'};

period_range = [10, 86400 + 3600];

dock on;figure(1);close all;

hao = 0;
if hao
    figure();
        topts = copts;
        topts.type = 'original';
        tsplot(TFs{1},topts);
        colororder(gcf,[0,0,0;0.5,0.5,0.5]);

    figure();
        topts = copts;
        topts.type  = 'error';
        topts.time_range = {'2003-10-29T06:00:00.000',...
                            '2003-10-29T09:00:00.000'};
        %tsplot(TFs,topts);

    figure();
        zopts = copts;
        zopts.type = 2;
        zplot(TFs,zopts);
    return
end

%% Time series plots
figure();
    topts = copts;
    topts.type = 'original';
    tsplot(TFs{1},topts);

if 0
    topts = copts;
    topts.type  = 'error';
    %topts.time_range =
    for i = 1:length(TFs)
        figure();
        tsplot(TFs{i},topts);
    end
end

if 1
    figure();
        topts = copts;
        topts.type  = 'error';
        tsplot(TFs,topts);
end

if (0)
    %% Compare
    figure();
        topts = copts;
        topts.printname = 'ts-error-tf1-tf2-tf3';
        topts.type  = 'error';
        topts.print = 1;
        tsplot(TFs,topts);


    figure();
        topts = copts;
        topts.printname = 'ts-error-zoom-tf1-tf2-tf3';
        topts.type  = 'error';
        topts.time_range = {'2003-10-29T06:00:00.000',...
                            '2003-10-29T09:00:00.000'};
        tsplot(TFs,topts);
end

% Compare
figure();
    snopts = copts;
    snopts.printname = 'sn-tf1-tf2-tf3';
    snopts.period_range = period_range;
    snplot(TFs,snopts);

%% Z plots
% Compare
figure();
    zopts = copts;
    zopts.type = 1;
    zopts.printname = 'z-tf1-tf2-tf3';
    zopts.period_range = period_range;
    zplot(TFs,zopts);

%% DFT plots
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftopts.period_range = period_range;
    dftplot(TFs{1},dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-phases';
    dftopts.period_range = period_range;
    dftplot(TFs{1},dftopts);

%%
figure();
    dftopts = copts;
    dftopts.type = 'error-averaged';
    dftopts.period_range = period_range;
    dftplot({TFs{1},TFs{2}},dftopts);

if 0
    figure()
        qqplot_(TF1,{},10,1,1);
end

if print_figs == 1
    figHTML(copts.printOptions.printDir)
end
