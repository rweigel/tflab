addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

chainid = 'KAP03';

siteid = 'KAP103';
daterange = '20031108-20031205';

siteid = 'KAP163';
daterange = '20031028-20031124';
daterange = '20031028-20031102';

if ~exist('TF1','var')
    fname = fullfile(scriptdir(),'data',chainid,siteid,...
                sprintf('%s-%s-tf1.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    TF1 = load(fname);
    TF1 = tflab_preprocess(TF1,'both',0,1);
    TF1 = tflab_metrics(TF1);    
end
if ~exist('TF2','var')
    fname = fullfile(scriptdir(),'data',chainid,siteid,...
                sprintf('%s-%s-tf2.mat',siteid, daterange));
    logmsg('Reading %s\n',fname);
    TF2 = load(fname);
    TF2 = tflab_preprocess(TF2,'both',0,1);
    TF2 = tflab_metrics(TF2);
end

%% Common plot options
copts.printdir = fullfile(scriptdir(),'data',siteid,'figures');
copts.print = 0;  % Set to 1 to print pdf of each figure created.
copts.printfmt = {'pdf'};

dock on;figure(1);close all;

%% Time series plots
figure();
    topts = copts;
    topts.type = 'original';
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

    
%% Z plots
% Compare
figure();
    zopts = copts;
    zopts.period_range = [1, 2*1e3];
    zopts.type = 2;
    zplot({TF1,TF2},zopts);

if 0
    
%% SN plots
figure();
    snopts = copts;
    snplot(TF1,snopts);

%%
figure();
    snopts = copts;
    snplot(TF2,snopts);
 
%%
% Compare
figure();
    snopts = copts;
    snplot({TF1,TF2},snopts);


%% DFT plots
figure();
    dftopts = copts;
    dftplot(TF1,dftopts);

%%
figure();
    dftopts = copts;
    dftopts.type = 'error-averaged';
    dftplot({TF1,TF2},dftopts);


%% Z plots
% Compare
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({TF1,TF2},zopts);
    
figure()
    qqplot_(TF1,{},10,1,1);
end