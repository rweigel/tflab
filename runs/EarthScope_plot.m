clear;
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';

start = '20160610';
stop = '20160616';

%start = '20160610';
%stop = '20160618';

%start = '20160610';
%stop = '20160623';

%start = '20160617';
%stop = '20160620';

outdir = fullfile(scriptdir(),'data','EarthScope',id);

for tfn = 1:3
    f{tfn} = fullfile(outdir, sprintf('VAQ58-%s-%s-tf%d.mat',start,stop,tfn));
end
TF1 = loadtf(f{1});
TF1 = tflab_preprocess(TF1);
TF1 = tflab_metrics(TF1);

TF2 = loadtf(f{2});
TF2 = tflab_preprocess(TF2);
TF2 = tflab_metrics(TF2);

TF3 = loadtf(f{3});
TF3 = tflab_preprocess(TF3);
TF3 = tflab_metrics(TF3);

%% Set common print options
copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(outdir,'figures');
copts.printfmt = {'pdf'};

dock on;figure(1);close all;

%% Time series plots

% Plot original time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsplot(TF1,tsopts);

if 0    
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

% Plot error for TF3 only
figure();
    tsopts = copts;
    tsopts.type = 'error';
    tsplot(TF3,tsopts);
end

% Compare all errors
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsplot({TF1,TF2, TF3},tsopts);

%% DFT plots
% Plot DFTs for TF1 only (will be same for both)
if 0
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TF1,dftopts);


figure();
    dftopts = copts;
    dftopts.type = 'error-averaged-magphase';
    dftplot(TF1,dftopts);
end


figure()
    y = abs(TF2.Out_.Predicted(:,1));
    [N,X] = hist(y(y < 500),20);
    semilogy(X,N/sum(N),'.','MarkerSize',20);
    hold on;grid on;
    y = abs(TF3.Out_.Predicted(:,1));
    [N3,X3] = hist(y(y < 500),20);
    semilogy(X3,N3/sum(N3),'.','MarkerSize',20);
    ylabel('Probability')
    xlabel('$E_x$ [mV/m]','Interpreter','Latex')

%% SN plots
figure();
    snopts = copts;
    %snopts.print = 1;
    snplot(TF1,snopts)

    
% Compare all
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snplot({TF1,TF3},snopts);

%% Z plots
if 0
figure();
    zopts = copts;
    zopts.type = 2;
    zplot(TF1,zopts);

figure();    
    zopts = copts;
    zopts.type = 2;
    zplot(TF3,zopts);
end

if 1
% Compare Z between TF1 and TF2    
figure();
    zopts = copts;
    zopts.period_range = [6,6*3600];
    zplot({TF1,TF3},zopts);
end

%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qopts = copts;
    qqplot_({TF1,TF2, TF3},qopts,fidx,comp,sidx);
