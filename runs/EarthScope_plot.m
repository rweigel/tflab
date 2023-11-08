clear;
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';

start = '20160610';
stop = '20160616';
%stop = '20160618';

outdir = fullfile(scriptdir(),'data','EarthScope',id);

for tfn = 1:3
    f{tfn} = fullfile(outdir, sprintf('VAQ58-%s-%s-tf%d.mat',start,stop,tfn));
end
TF1 = loadtf(f{1});
TF1 = tflab_preprocess(TF1);
TF1 = tflab_metrics(TF1);

TF2 = loadtf(f{2});
TF2 = tflab_preprocess(TF2);
TF2 = tflab_metrics(TF2);

TF3 = loadtf(f{3});
TF3 = tflab_preprocess(TF3);
TF3 = tflab_metrics(TF3);

%% Set common print options
copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};

figure(1);close all;

%% Time series plots

% Plot original time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsplot(TF1,tsopts);

if 0    
% Plot error for TF1 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(TF1,tsopts);

% Plot error for TF2 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(TF2,tsopts);

% Plot error for TF3 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(TF3,tsopts);
end

% Compare all errors
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({TF1,TF2, TF3},tsopts);

%% DFT plots
% Plot DFTs for TF1 only (will be same for both)
if 0
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TF1,dftopts);


figure();
    dftopts = copts;
    dftopts.type = 'error-averaged-magphase';
    dftplot(TF1,dftopts);
end

%% SN plots
% Compare SN between TF1 and TF2
figure();
    snopts = copts;
    snplot({TF1,TF2,TF3},snopts);

%% Z plots
if 0
figure();
    zopts = copts;
    zopts.type = 2;
    zplot(TF1,zopts);

figure();    
    zopts = copts;
    zopts.type = 2;
    zplot(TF3,zopts);
end

if 1
% Compare Z between TF1 and TF2    
figure();
    zopts = copts;
    zopts.period_range = [3, 86400];
    zplot({TF1,TF2,TF3},zopts);
end

%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qopts = copts;
    qqplot_({TF1,TF2, TF3},qopts,fidx,comp,sidx);
