clear;
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

% Set common print options

short_run = 0;
print_figs = 0;

if short_run
    start = '20120712';
    stop = '20120716';
else
    start = '20120712';
    stop = '20121107';
end

dirstr  = sprintf('tfs-%s-%s',start,stop);
rundir = fullfile(scriptdir(),'data','Middelpos',dirstr);

for tfn = 1:4
    fname{tfn} = fullfile(rundir, sprintf('Middelpos-tf%d.mat',tfn));
    TFs{tfn} = loadtf(fname{tfn});
    %TFs{tfn} = tflab_preprocess(TFs{tfn});
    %TFs{tfn} = tflab_metrics(TFs{tfn});
end

copts = struct();
    copts.print = print_figs; % Set to 1 to print pdf of each figure created.
    copts.printOptions.printDir = fullfile(rundir,'figures');
    copts.printOptions.printFormats = {'pdf'};

dock on;figure(1);close all;

%% Time series plots
tsopts = copts;
tsopts.printOptions.printFormats = {'pdf'};

if (1)
    figure();
        tsopts.type = 'original';
        tsplot(TFs{1},tsopts);
    figure();
        tsopts.type = 'final';
        tsplot(TFs{3},tsopts);
end

if (1)
    tsopts.type = 'error';
    figure();
        tsplot(TFs{1},tsopts);
    figure();
        tsplot(TFs{3},tsopts);
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

if 0
    figure();
        dftopts = copts;
        dftopts.type = 'original-averaged-reals';
        dftplot(TFs{1},dftopts);

    figure();
        dftopts = copts;
        dftopts.type = 'original-averaged-imaginaries';
        dftplot(TFs{1},dftopts);

    % Plot Fourier phases for In/Out of TF1 (will be the same for both)
    figure();
        dftopts = copts;
        dftopts.type = 'original-averaged-phases';
        dftplot(TFs{1},dftopts);

        % Compare TF1 and TF2
        figure();
            dftopts = copts;
            dftopts.type = 'error-averaged-realimag';
            dftplot(TFs{1},dftopts);
end

%% SN plots
figure();
    snopts = copts;
    snplot(TFs(:),snopts);

%% Z plots
figure();
    zopts = copts;
    zopts.type = 1;
    %zopts.period_range = [1, 86400];
    zplot(TFs,zopts);

figure();
    qqopts = copts;
    qqopts.printOptions.printDir = fullfile(rundir,'figures','qqplot');
    fidx = 20; % frequency number
    comp = 2;  % component (x = 1, y = 2)
    qqplot_(TFs{1},qqopts,comp,fidx);

figure();
        qqopts = copts;
        qqopts.printOptions.printDir = fullfile(rundir,'figures','qqplot');
        qqopts.type = 'combined';
        qqplot_(TFs{1},qqopts);
