
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

%% Load data if not already in memory.
daterange = '20120712-20120721';
%daterange = '20120712-20121031';
if ~exist('TF1','var')
    fname = fullfile('data','Middelpos',...
            sprintf('Middelpos-%s-tf1.mat',daterange));
    fnamefull = fullfile(scriptdir(),fname);
    logmsg(sprintf('Reading %s',fname));
    TF1 = load(fnamefull);
    TF1 = tflab_preprocess(TF1,'both',0);
    TF1 = tflab_metrics(TF1);
end
if ~exist('TF2','var')
    fname = fullfile('data','Middelpos',...
            sprintf('Middelpos-%s-tf2.mat',daterange));
    fnamefull = fullfile(scriptdir(),fname);
    logmsg(sprintf('Reading %s',fname));
    TF2 = load(fnamefull);
    TF2 = tflab_preprocess(TF2,'both',0);
    TF2 = tflab_metrics(TF2);
end

%% Set common print optionsx

copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};

dock on;figure(1);close all;

%% Time series plots
% Plot raw time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsplot(TF1,tsopts);

% Plot error for S1 only
if 0
    figure();
        tsopts = copts;
        tsopts.type = 'error';
        tsplot(TF1,tsopts);
end

 % Compare S1 and S2 error
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({TF1,TF2},tsopts);


%% DFT plots
% Plot Fourier amplitudes for In/Out of TF1 (will be the same for both)
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';    
    dftplot(TF1,dftopts);

if 0    
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-reals';
    dftplot(TF1,dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-imaginaries';
    dftplot(TF1,dftopts);
end

% Plot Fourier phases for In/Out of TF1 (will be the same for both)
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-phases';
    dftplot(TF1,dftopts);

% Compare TF1 and TF2
figure();
    dftopts = copts;
    dftopts.type = 'error-averaged';
    dftplot({TF1,TF2},dftopts);

    
%% SN plots
% Plot SN for TF1 only
figure();
    snopts = copts;
    snplot(TF1,snopts);

% Plot SN for TF2 only
figure();
    snopts = copts;
    snplot(TF2,snopts);

% Compare SN between TF1 and TF2
figure();
    snopts = copts;
    snplot({TF1,TF2},snopts);

%% Z plots
% Compare Z between S1 and S2
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({TF1,TF2},zopts);


%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 20; % frequency number
comp = 2;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qqplot_(TF2,fidx,comp,sidx);
    