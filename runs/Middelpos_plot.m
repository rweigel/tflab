%%
addpath(fullfile(scriptdir(),'..','plot'));

daterange = '20120712-20121031';
if ~exist('S1','var')
    load(fullfile(scriptdir(),'data','Middelpos',...
        sprintf('Middelpos-%s-tf1.mat',daterange)));
end
if ~exist('S2','var')
    load(fullfile(scriptdir(),'data','Middelpos',...
        sprintf('Middelpos-%s-tf2.mat',daterange)));
end
if ~exist('S3','var')
    %load(fullfile(scriptdir(),'data','Middelpos',sprintf('Middelpos-20120712-20121031-tf3.mat',datarange)));
end

%% Common options
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
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
    