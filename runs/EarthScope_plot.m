clear;
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

id = 'VAQ58';
%id = 'ORF03';
outdir = fullfile(scriptdir(),'data','EarthScope',id);

%% Set common print options
copts.print    = 0; % Set to 1 to print
copts.printdir = fullfile(outdir,'figures');
copts.printfmt = {'pdf','png'};

if 0
    start = '20160610';
    stop = '20160614';
    time_range_full = {'2016-06-11T00:00:00.000','2016-06-16T12:00:00.000'};
    time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-15T06:00:00.000'};
end

if 1
    start = '20160610';
    stop = '20160616';
    time_range_full = {'2016-06-11T00:00:00.000','2016-06-16T12:00:00.000'};
    time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-15T06:00:00.000'};
end

if 0
    start = '20160610';
    stop = '20160623';
    time_range_full = {'2016-06-11T00:00:00.000','2016-06-16T12:00:00.000'};
    time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-23T12:00:00.000'};
end

if strcmp(id,'ORF03')
    start = '20070831';
    stop = '20070903';
    time_range_full = {'2007-08-31T00:00:00.000','2007-09-04T00:00:00.000'};
    time_range_zoom = time_range_full;
end

for tfn = 1:3
    fname = fullfile(outdir, sprintf('%s-%s-%s-tf%d.mat',id,start,stop,tfn));
    TFs{tfn} = loadtf(fname,'',1,1);
end

dock on;figure(1);close all;

%% Time series plots

% Plot original time series data used for TF1 (will be same as that for TF2)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsopts.printname = 'ts-tf1';
    tsopts.print = 0;
    tsplot(TFs{1},tsopts);

if 1
    % Plot error for TF1 only
    figure();
        tsopts = copts;
        tsopts.type = 'error';
        tsplot(TFs{1},tsopts);
    
    % Plot error for TF2 only
    figure();
        tsopts = copts;
        tsopts.type = 'error';
        tsplot(TFs{2},tsopts);
    
    % Plot error for TF3 only
    figure();
        tsopts = copts;
        tsopts.type = 'error';
        tsplot(TFs{3},tsopts);
end

%%
% Compare errors
figure();
    tsopts = copts;
    tsopts.time_range = time_range_full;
    tsopts.type  = 'error';
    tsopts.printname = 'ts-error-tf1-tf3';
    tsplot({TFs{1},TFs{3}},tsopts);

figure();
    tsopts = copts;
    tsopts.time_range = time_range_zoom;
    tsopts.type  = 'error';
    tsopts.printname = 'ts-error-zoom-tf1-tf3';
    tsplot({TFs{1},TFs{3}},tsopts);


%% DFT plots
% Plot DFTs for TF1 only (will be same for both)
if 0
    figure();
        dftopts = copts;
        dftopts.type = 'original-averaged';
        dftplot(TFs{1},dftopts);
    
    
    figure();
        dftopts = copts;
        dftopts.type = 'error-averaged-magphase';
        dftplot(TFs{2},dftopts);
end

%%
figure()

    fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
    time_range = {'2016-06-14T18:00:00.000','2016-06-15T06:00:00.000'};
    mldn_range(1) = datenum(time_range{1},fmt);
    mldn_range(2) = datenum(time_range{2},fmt);
    to = datenum(TFs{2}.Metadata.timestart,fmt);
    ppd = 86400;
    nt = size(TFs{2}.Out_.Predicted,1);
    t = to + (0:nt-1)'/ppd;
    tidx = find(t >= mldn_range(1) & t <= mldn_range(2));

    figprep();    
    popts = tflabplot_options(TFs{2}, copts, 'tsplot');
    popts.Positions = {popts.PositionTop, popts.PositionBottom};
    popts.Positions{1}(4) = 0.35;
    popts.Positions{2}(4) = 0.35;
    popts.Positions{1}(2) = 0.49;
    for comp = 1:2
    %ax(comp) = subplot('Position',popts.Positions{comp});
    subplot(2,1,comp)
        y = abs(TFs{2}.Out_.Predicted(tidx,1));
        [N,X] = hist(y,20);
        semilogy(X,N/sum(N),'.','MarkerSize',20);
        hold on;grid on;
        y = abs(TFs{3}.Out_.Predicted(tidx,comp));
        [N3,X3] = hist(y,20);
        semilogy(X3,N3/sum(N3),'.','MarkerSize',20);
        legend('TFLab','EMTF')
        ylabel('Probability')
        xlabel(sprintf('$%s$ Predicted [mV/m]',TFs{2}.Metadata.outstr{comp}),'Interpreter','Latex')
    end
    if copts.print == 1   
        figsave(fullfile(copts.printdir, 'pdf-tf1-tf3.png'));
    end

%% SN plots
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snopts.printname = 'sn-tf1';
    snopts.print = 0;
    snplot(TFs{1},snopts,1);

% Compare all
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snopts.printname = 'sn-tf1-tf3';
    snopts.print = 0;
    snplot(TFs,snopts);

% Compare
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snopts.printname = 'sn-tf1-tf3';
    snopts.print = 0;
    snplot({TFs{1},TFs{3}},snopts);

%% Z plots
if 0
    figure();
        zopts = copts;
        zopts.type = 2;
        zplot(TFs{1},zopts);
    
    figure();    
        zopts = copts;
        zopts.type = 2;
        zplot(TFs{2},zopts);
end

if 1
    % Compare Z
    figure();
        zopts = copts;
        %zopts.print = 1;
        zopts.printname = 'z-tf1-tf3';
        zopts.period_range = [7,6*3600];
        zopts.unwrap = 0;
        zopts.type = 1;
        zplot({TFs{2},TFs{3}},zopts);
        %zplot({TFs{1},TFs{2},TFs{3}},zopts);
end

if 0
    % Should match http://ds.iris.edu/spud/emtf/15014571
    zopts.type = 2;
    figure();
        zplot(TFs,zopts,2);
    figure();
        zplot(TFs,zopts,3);
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
    %qqplot_({TF1,TF2,TF3},qopts,fidx,comp,sidx);
    qqplot_({TFs{1},TFs{3}},qopts,fidx,comp,sidx);
