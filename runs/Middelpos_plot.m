
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

%% Load data if not already in memory.
%daterange = '20120712-20120717';
daterange = '20120712-20121107';
if ~exist('TF1','var')
    fname = fullfile('data','Middelpos',...
            sprintf('Middelpos-%s-tf1.mat',daterange));
    fnamefull = fullfile(scriptdir(),fname);
    logmsg(sprintf('Reading %s',fname));
    TF1 = load(fnamefull);
end
if ~exist('TF2','var')
    fname = fullfile('data','Middelpos',...
            sprintf('Middelpos-%s-tf2.mat',daterange));
    fnamefull = fullfile(scriptdir(),fname);
    logmsg(sprintf('Reading %s',fname));
    TF2 = load(fnamefull);
end

%% Set common print optionsx

copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};

figure(1); 
close all;

%% Time series plots
% Plot raw time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'raw';
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

%% SN plots
if 0
    % Plot SN for S1 only
    figure();
        snopts = copts;
        %snopts.period_range = [1, 86400];
        snplot(TF1,snopts);

    % Plot SN for S2 only
    figure();
        snopts = copts;
        %snopts.period_range = [1, 86400];
        snplot(TF2,snopts);
end

% Compare SN between S1 and S2
figure();
    snopts = copts;
    %snopts.period_range = [1, 110*86400];
    snplot({TF1,TF2},snopts);


%% PSD plots
% Plot PSDs for S1 only (will be the same for both)
figure();
    psdopts = copts;
    %psdopts.period_range = [1, 86400];
    psdopts.type = 'smoothed';    
    psdplot(TF1,psdopts);

% Compare PSD between S1 and S2
figure();
    psdopts = copts;
    psdopts.type = 'error-smoothed';
    %psdopts.period_range = [1, 86400];
    psdplot({TF1,TF2},psdopts);


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

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qqplot_(TF1,fidx,comp,sidx);
    