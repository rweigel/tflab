%%
addpath(fullfile(scriptdir(),'..','plot'));


%% Common options
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.print = 1;
copts.printfmt = {'pdf'};


%%
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



%%
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
    snplot({S1,S2},sopts);


%%
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot(S1,psdopts);

%%
figure();
    psdopts = copts;
    psdopts.period_range = [1, 86400];
    psdplot({S1,S2},psdopts);


%%
% Compare
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot({S1,S2},zopts);