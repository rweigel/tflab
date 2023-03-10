addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';

%% Load data if not already in memory.
daterange = '20160610-20160624';
%daterange = '20120712-20120717';
if ~exist('TF1','var')
    fname = fullfile('data','EarthScope',id,...
            sprintf('%s-%s-tf1.mat',id,daterange));
    fnamefull = fullfile(scriptdir(),fname);
    logmsg(sprintf('Reading %s',fname));
    TF1 = load(fnamefull);
end

%% Set common print options
copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};


%% Time series plots
% Plot raw time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'raw';
    tsplot(TF1,tsopts);

% Plot error for S1 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(TF1,tsopts);

if 0    
% Compare S1 and S2 error
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({TF1,TF2},tsopts);
end    

%% SN plots
% Plot SN for S1 only
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(TF1,snopts);

if 0    
% Plot SN for S2 only
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(TF2,snopts);
 
% Compare SN between S1 and S2
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot({TF1,TF2},snopts);
end

%% PSD plots
% Plot PSDs for S1 only
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot(TF1,psdopts);

% 
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot(TF1,zopts);
    
if 0    
% Compare SN between S1 and S2
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot({TF1,TF2},psdopts);


%% Z plots
% Compare Z between S1 and S2
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
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
    