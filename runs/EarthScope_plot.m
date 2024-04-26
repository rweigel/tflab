function EarthScope_plot(id, print_figs)

if ~exist('print_figs','var') || print_figs ~= 1
    print_figs = 0;
end

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

if strcmp(id,'VAQ58')
    start = '2016-06-10T18:19:12';
    stop = '2016-06-14T18:19:11';
    time_range_full = {start, stop};
    time_range_zoom = {};
    if 0
        start = '20160610';
        stop = '20160614';
        time_range_full = {start, stop};
        time_range_zoom = {};
    end

    if 0
        start = '20160610';
        stop = '20160616';
        time_range_full = {'2016-06-11T00:00:00.000','2016-06-16T12:00:00.000'};
        time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-15T06:00:00.000'};
    end

    if 0
        start = '20160610';
        stop = '20160623';
        time_range_full = {'2016-06-11T00:00:00.000','2016-06-16T12:00:00.000'};
        time_range_zoom = {};
        %time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-23T12:00:00.000'};
    end
end

if strcmp(id,'ORF03') || strcmp(id,'ORG03')
    %start = '2007-08-19T01:48:36';
    %stop = '2007-09-07T17:18:40';
    start = '2007-08-31T01:48:36';
    stop = '2007-09-04T01:48:35';
    time_range_full = {start, stop};
    time_range_zoom = {};
end

dno = datenum(start,'yyyy-mm-ddTHH:MM:SS');
dnf = datenum(stop,'yyyy-mm-ddTHH:MM:SS');
dso = datestr(dno,'yyyymmddTHHMMSS');
dsf = datestr(dnf,'yyyymmddTHHMMSS');
dirstr  = sprintf('tfs-%s-%s',dso,dsf);
rundir = fullfile(scriptdir(),'data','EarthScope',id,dirstr);

%% Set common print options
copts.print = print_figs; % Set to 1 to print pdf of each figure created.
copts.printOptions.printDir = fullfile(rundir,'figures');
copts.printOptions.printFormats = {'png'};

for tfn = 1:3
    fname = fullfile(rundir, sprintf('%s-tf%d.mat',id,tfn));
    TFs{tfn} = loadtf(fname);
end

if print_figs
    dock off;close all;
else
    dock on;figure(1);close all;
end

%% Time series plots

% Plot original time series data used for TF1 (will be same for all)
figure();
    tsopts = copts;
    tsopts.type = 'original';
    tsplot(TFs{1},tsopts);

% Compare errors
figure();
    tsopts = copts;
    tsopts.type  = 'error';
    tsopts.printOptions.printName = 'tsplot-error-tf1-tf3';
    tsplot({TFs{1},TFs{3}},tsopts);
    if ~isempty(time_range_zoom)
        figure();
            tsopts = copts;
            tsopts.time_range = time_range_zoom;
            tsopts.type  = 'error';
            tsopts.printname = 'tsplot-error-zoom-tf1-tf3';
            tsplot({TFs{1},TFs{3}},tsopts);
    end

%% Z plots
if 1
    figure();
        zopts = copts;
        zopts.type = 1;
        zplot(TFs{1},zopts);

    figure();
        zopts = copts;
        zopts.type = 1;
        zplot(TFs{3},zopts);
end

if 1
    % Compare
    figure();
        zopts = copts;
        %zopts.print = 1;
        zopts.printname = 'z-tf1-tf3';
        zopts.period_range = [7,6*3600];
        zopts.unwrap = 0;
        zopts.type = 1;
        %zplot(TFs,zopts);
        zplot({TFs{1},TFs{3}},zopts);
end

%% SN plots
% Compare all
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snopts.printname = 'sn-tf1-tf3';
    snplot({TFs{1},TFs{3}},snopts);


%% DFTs
if 1
    dftopts = copts;
    dftopts.period_range = [7,24*3600];
    figure();
        dftopts.type = 'original-averaged';
        dftplot(TFs{1},dftopts);

    figure();
        dftopts = copts;
        dftopts.type = 'original-averaged-phases';
        dftplot(TFs{1},dftopts);

    figure();
        dftopts.type = 'error-averaged-magphase';
        dftplot(TFs{1},dftopts);
end

%% Histograms
if 0
    figure();

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
end

if 0
    % Should match
    % VAQ58: http://ds.iris.edu/spud/emtf/15014571
    % ORF03: http://ds.iris.edu/spud/emtf/14866915
    % When adding 180 degrees
    % But ... if phase is negative, they plot as postive ... or not at all
    % https://ds.iris.edu/spudservice/data/15014347
    zopts.type = 2;
    figure();
        zplot(TFs{3},zopts,2);
    figure();
        zplot(TFs{3},zopts,3);
end

fidx = 10; % frequency index
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment index

figure();
    qopts = copts;
    qqplot_(TFs{3},qopts,comp,fidx);
figure();
    qopts = copts;
    qqplot_(TFs{1},qopts,comp,fidx);

if print_figs == 1
    figHTML(copts.printOptions.printDir)
end


