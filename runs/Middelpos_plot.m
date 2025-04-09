function Middelpos_plot(rundir, filestr, print_figs)

logmsg('\n');
for tfn = 1:5
    fname{tfn} = fullfile(rundir, sprintf('%s-tf%d.mat',filestr,tfn));
    TFs{tfn} = loadtf(fname{tfn});
end

outdir = fullfile(scriptdir(),'data','KAP03','KAP103','tfs');
f = fullfile(outdir, sprintf('%s-%s.mat','KAP103','unknown-unknown-tf3'));
TFs{6} = loadtf(f);
%keyboard
TFs{6}.In = TFs{1}.In;
TFs{6}.Out = TFs{1}.Out;
TFs{6}.Z = -TFs{6}.Z;
TFs{6} = tflab_preprocess(TFs{6});
TFs{6} = tflab_metrics(TFs{6});

f = fullfile(outdir, sprintf('%s-%s.mat','KAP103','20031108-20031118-tf1'));
TFs{7} = loadtf(f);
TFs{7}.In = TFs{1}.In;
TFs{7}.Out = TFs{1}.Out;
TFs{7}.Z = -TFs{7}.Z;
TFs{7} = tflab_preprocess(TFs{7});
TFs{7} = tflab_metrics(TFs{7});

copts = struct();
    copts.print = print_figs; % Set to 1 to print pdf of each figure created.
    copts.printOptions.printDir = fullfile(rundir,'figures');
    copts.printOptions.printFormats = {'pdf','png'};
    copts.title = '';

keyboard
    
figure();
    zopts = copts;
    zopts.type = 1;
    %zopts.period_range = [1, 86400];
    zplot({TFs{1}, TFs{6}, TFs{7}},zopts);


figure();
    snopts = copts;
    snplot(TFs,snopts);

if print_figs
    dock off;close all;
    if exist(copts.printOptions.printDir, 'dir')
        rmdir(copts.printOptions.printDir,'s');
    end
    mkdir(copts.printOptions.printDir);
else
    dock on;figure(1);close all;
end

%% Time series plots
tsopts = copts;
if (1)
    figure();
        tsopts.type = 'original';
        tsplot(TFs{1},tsopts);
    figure();
        tsopts.type = 'detrended';
        tsplot(TFs{5},tsopts);
end

if (1)
    tsopts = copts;
    tsopts.type = 'error';
    for i = 1:length(TFs)
        figure();
            tsplot(TFs{i},tsopts);
    end
end

%% DFT plots
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TFs{1},dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-phases';
    dftplot(TFs{1},dftopts);

%figure();
    %dftopts.type = 'error-averaged-magphase';
    %dftplot(TFs,dftopts);

%% SN plots
figure();
    snopts = copts;
    snplot(TFs([1,3]),snopts);

figure();
    snopts = copts;
    snplot(TFs,snopts);
    
%% Z plots
figure();
    zopts = copts;
    zopts.type = 1;
    %zopts.period_range = [1, 86400];
    zplot(TFs([1,3]),zopts);

figure();
    zopts = copts;
    zopts.type = 1;
    %zopts.period_range = [1, 86400];
    zplot(TFs,zopts);


%% qq plots    
figure();
    qqopts = copts;
    %qqopts.printOptions.printDir = fullfile(rundir,'figures','qqplot');
    fidx = 20; % frequency number
    comp = 2;  % component (x = 1, y = 2)
    qqplot_(TFs{1},qqopts,comp,fidx);

figure();
    qqopts = copts;
    qqopts.type = 'combined';
    %qqplot_(TFs,qqopts);

if print_figs == 1
    figHTML(copts.printOptions.printDir)
end
