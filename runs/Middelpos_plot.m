
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

%% Load data if not already in memory.
daterange = '20120712-20121031';
%daterange = '20120712-20120717';
if ~exist('S1','var')
    fname = fullfile(scriptdir(),'data','Middelpos',...
                sprintf('Middelpos-%s-tf1.mat',daterange));
    logmsg(sprintf('Reading %s',fname));
    S1 = load(fname);
end
if ~exist('S2','var')
    fname = fullfile(scriptdir(),'data','Middelpos',...
                sprintf('Middelpos-%s-tf2.mat',daterange));
    logmsg(sprintf('Reading %s',fname));
    S2 = load(fname);
end

%% Set common print options
copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};


%% Time series plots
% Plot raw time series data used for S1 (will be same as that for S2)
figure();
    tsopts = copts;
    tsopts.type = 'raw';
    tsplot(S1,tsopts);

% Plot error for S1 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(S1,tsopts);

% Compare S1 and S2 error
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({S1,S2},tsopts);



%% SN plots
% Plot SN for S1 only
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(S1,snopts);

% Plot SN for S2 only
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(S2,snopts);
 
% Compare SN between S1 and S2
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot({S1,S2},snopts);


%% PSD plots
% Plot PSDs for S1 only
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot(S1,psdopts);

% Compare SN between S1 and S2
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot({S1,S2},psdopts);


%% Z plots
% Compare Z between S1 and S2
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({S1,S2},zopts);


%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qqplot_(S1,fidx,comp,sidx);
    