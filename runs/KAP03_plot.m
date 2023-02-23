%%
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

siteid = 'KAP103';
daterange = '20031108-20031205';
if ~exist('S1','var')
    fname = fullfile(scriptdir(),'data',siteid,...
                sprintf('%s-%s-tf1.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    S1 = load(fname);
end
if ~exist('S2','var')
    fname = fullfile(scriptdir(),'data',siteid,...
                sprintf('%s-%s-tf1.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    S2 = load(fname);
end

%% Common plot options
copts.printdir = fullfile(scriptdir(),'data',siteid,'figures');
copts.print = 0;
copts.printfmt = {'pdf'};

%% Time series plots
figure();
    topts = copts;
    topts.type = 'raw';
    tsplot(S1,topts);

%%
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot(S1,topts);

%% Compare
figure();
    topts = copts;
    topts.type  = 'error';
    tsplot({S1,S2},topts);


%% SN plots
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(S1,snopts);

%%
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot(S2,snopts);
 
%%
% Compare
figure();
    snopts = copts;
    snopts.period_range = [1, 86400];
    snplot({S1,S2},snopts);


%% PSD plots
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot(S1,psdopts);

%%
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot({S1,S2},psdopts);


%% Z plots
% Compare
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({S1,S2},zopts);
    
figure()
    qqplot_(S1,10,1,1);