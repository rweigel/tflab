addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';


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
    
% Compare TF1 and TF2 error
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({TF1,TF2},tsopts);

%% DFT plots
% Plot DFTs for TF1 only (will be same for both)
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TF1,dftopts);

%% SN plots
% Compare SN between TF1 and TF2
figure();
    snopts = copts;
    snplot({TF1,TF2},snopts);

%% Z plots
figure();
    zopts = copts;
    zopts.type = 1;
    zplot(TF1,zopts);

if 0
% Compare Z between TF1 and TF2    
figure();
    zopts = copts;
    zopts.period_range = [3, 86400];
    zplot({TF1,TF2},zopts);
end

%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qqplot_(TF1,fidx,comp,sidx);
    