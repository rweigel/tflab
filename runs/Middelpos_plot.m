
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();
outdir = fullfile(scriptdir(),'data','Middelpos');

start = '20120712';
stop = '20120721';

for tfn = 1:1
    f{tfn} = fullfile(outdir, sprintf('Middelpos-%s-%s-tf%d.mat',start,stop,tfn));
    TFs{tfn} = loadtf(f{tfn});
    TFs{tfn} = tflab_preprocess(TFs{tfn});
    TFs{tfn} = tflab_metrics(TFs{tfn});
end

%% Set common print options

copts.print    = 0; % Set to 1 to print pdf of each figure created.
copts.printdir = fullfile(scriptdir(),'data','Middelpos','figures');
copts.printfmt = {'pdf'};

dock on;figure(1);close all;

%% Time series plots
% Plot raw time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsplot(TFs{1},tsopts);

    
% Plot error for S1 only
if 1
    figure();
        tsopts = copts;
        tsopts.type = 'error';
        tsplot(TFs{1},tsopts);
end

if 0
    % Compare S1 and S2 error
    figure();
        tsopts = copts;
        tsopts.type  = 'error';
        %tsplot({TF1,TF2},tsopts);
end

%% DFT plots
% Plot Fourier amplitudes for In/Out of TF1 (will be the same for both)
if 0    
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';    
    dftplot(TFs{1},dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-reals';
    dftplot(TFs{1},dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-imaginaries';
    dftplot(TFs{1},dftopts);

% Plot Fourier phases for In/Out of TF1 (will be the same for both)
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-phases';
    dftplot(TF{1},dftopts);

    % Compare TF1 and TF2
    figure();
        dftopts = copts;
        dftopts.type = 'error-averaged';
        dftplot(TF(1:2),dftopts);
end

    
%% SN plots
% Plot SN for TF1 only
if 0
    figure();
        snopts = copts;
        snplot(TF1,snopts);
    
    % Plot SN for TF2 only
    figure();
        snopts = copts;
        snplot(TF2,snopts);
end

% Compare all
figure();
    snopts = copts;
    snplot(TFs{1},snopts);

%% Z plots
% Compare Z between S1 and S2
figure();
    zopts = copts;
    zopts.period_range = [1, 86400];
    zplot(TFs(:),zopts);


%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 20; % frequency number
comp = 2;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qqplot_(TFs{1},struct(),fidx,comp,sidx);
    