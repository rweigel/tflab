addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

siteid = 'KAP103';
chainid = 'KAP03';

daterange = '20031108-20031205';
if ~exist('TF1','var')
    fname = fullfile(scriptdir(),'data',chainid,siteid,...
                sprintf('%s-%s-tf1.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    TF1 = load(fname);
end
if ~exist('S2','var')
    fname = fullfile(scriptdir(),'data',chainid,siteid,...
                sprintf('%s-%s-tf2.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    TF2 = load(fname);
end

%% Common plot options
copts.printdir = fullfile(scriptdir(),'data',siteid,'figures');
copts.print = 0;  % Set to 1 to print pdf of each figure created.
copts.printfmt = {'pdf'};

%% Time series plots
figure();
    topts = copts;
    topts.type = 'raw';
    tsplot(TF1,topts);

%%
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot(TF1,topts);

%% Compare
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot({TF1,TF2},topts);


%% SN plots
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(TF1,snopts);

%%
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(TF2,snopts);
 
%%
% Compare
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot({TF1,TF2},snopts);


%% PSD plots
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot(TF1,psdopts);

%%
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot({TF1,TF2},psdopts);


%% Z plots
% Compare
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({TF1,TF2},zopts);
    
figure()
    qqplot_(TF1,10,1,1);