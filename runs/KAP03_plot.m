addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

chainid = 'KAP03';

siteid = 'KAP103';
daterange = '20031108-20031205';

siteid = 'KAP163';
%daterange = '20031028-20031124';
%daterange = '20031028-20031102';
tfs = {...
        '20031028-20031031-tf1a',...
        '20031031-20031103-tf1b',...
        '20031028-20031124-tf2'};

outdir = fullfile(scriptdir(),'data','KAP03',siteid);

for tfn = 1:length(tfs)
    f = fullfile(outdir, sprintf('%s-%s.mat',siteid,tfs{tfn}));
    TFs{tfn} = loadtf(f);
    TFs{tfn} = tflab_preprocess(TFs{tfn});
    TFs{tfn} = tflab_metrics(TFs{tfn});
end

%% Common plot options
copts.printdir = fullfile(scriptdir(),'data',chainid,siteid,'figures');
copts.print = 0;  % Set to 1 to print pdf of each figure created.
copts.printfmt = {'pdf'};

dock on;figure(1);close all;

hao = 0;
if hao
    figure();
        topts = copts;
        topts.type = 'original';
        tsplot(TF1,topts);
        colororder(gcf,[0,0,0;0.5,0.5,0.5]);

    figure();
        topts = copts;
        topts.type  = 'error';
        topts.time_range = {'2003-10-29T06:00:00.000',...
                            '2003-10-29T09:00:00.000'};
        tsplot({TF1,TF2},topts);

    figure();
        zopts = copts;
        zopts.type = 2;
        zplot({TF1,TF2},zopts);
    return
end

%% Time series plots
figure();
    topts = copts;
    topts.type = 'original';
    tsplot(TFs{3},topts);

%%
if 0
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot(TFs{1},topts);


%% Compare
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot({TFs{1},TFs{2}},topts);


figure();
    topts = copts;
    topts.type  = 'error';
    topts.time_range = {'2003-10-29T06:00:00.000',...
                        '2003-10-29T09:00:00.000'};
    tsplot({TFs{1},TFs{2}},topts);
end

if 0
%% SN plots
figure();
    snopts = copts;
    snplot(TFs{1},snopts);

figure();
    snopts = copts;
    snplot(TFs{2},snopts);

% Compare
figure();
    snopts = copts;
    snplot({TFs{1},TFs{2}},snopts);
end

%% Z plots
if 1
% Compare
figure();
    zopts = copts;
    zopts.period_range = [10, 2*1e3];
    zopts.type = 1;
    zplot(TFs,zopts);
end
    
if 0
%% DFT plots
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TFs{1},dftopts);

%%
figure();
    dftopts = copts;
    dftopts.type = 'error-averaged';
    dftplot({TFs{1},TFs{2}},dftopts);
end

if 0    
figure()
    qqplot_(TF1,{},10,1,1);
end