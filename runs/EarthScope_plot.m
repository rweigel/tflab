clear;
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

print_figs = 1;

%id = 'VAQ58';
%id = 'ORG03';
id = 'ORF03';

if strcmp(id,'ORG03')
    start = '20070831';
    stop = '20070904';
    time_range_full = {'2007-08-31T00:00:00.000','2007-09-04T00:00:00.000'};
    time_range_zoom = {};
end

if strcmp(id,'ORF03')
    %start = '2007-08-19T01:48:36';
    %stop = '2007-09-07T17:18:40';
    start = '2007-08-31T01:48:36';
    stop = '2007-09-04T01:48:35';
    dno = datenum(start,'yyyy-mm-ddTHH:MM:SS');
    dnf = datenum(stop,'yyyy-mm-ddTHH:MM:SS');
    dirstr  = sprintf('tfs-%s-%s',...
        datestr(dno,'yyyymmddTHHMMSS'),datestr(dnf,'yyyymmddTHHMMSS'));
    rundir = fullfile(scriptdir(),'data','EarthScope',id,dirstr);
    time_range_full = {'2007-08-31T00:00:00.000','2007-09-04T00:00:00.000'};
    time_range_zoom = {};
end

if strcmp(id,'VAQ58')
    if 0
        start = '20160610';
        stop = '20160614';
        time_range_full = {'2016-06-11T00:00:00.000','2016-06-15T00:00:00.000'};
        time_range_zoom = {};
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
        time_range_zoom = {};
        %time_range_zoom = {'2016-06-14T18:00:00.000','2016-06-23T12:00:00.000'};
    end
end


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

if 1
    for tfn = 1:3
        % Plot original time series data used for TF1 (will be same for all)
        figure();
            tsopts = copts;
            tsopts.type = 'error';
            tsplot(TFs{tfn},tsopts);
    end
end

%% Compare errors
figure();
    tsopts = copts;
    %tsopts.time_range = time_range_full;
    tsopts.type  = 'error';
    tsopts.printOptions.printName = 'ts-error-tf1-tf3';
    tsplot({TFs{1},TFs{3}},tsopts);

    if ~isempty(time_range_zoom)
        figure();
            tsopts = copts;
            tsopts.time_range = time_range_zoom;
            tsopts.type  = 'error';
            tsopts.printname = 'ts-error-zoom-tf1-tf3';
            tsplot({TFs{1},TFs{3}},tsopts);
    end

%% DFTs
% Plot DFTs for TF1 only (will be same for both)
if 0
figure();
    dftopts = copts;
    dftopts.type = 'error-averaged';
    dftplot({TFs{1}, TFs{3}},dftopts);
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

%% SN plots
% Compare all
figure();
    snopts = copts;
    snopts.period_range = [7,6*3600];
    snopts.printname = 'sn-tf1-tf3';
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

if 1
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

%% Regression plots
% Plot regression errors for a component of S1's Z at a single frequency
% for one of the segments. (For S2, there is only one segment that was
% used to compute Z.)

fidx = 10; % frequency number
comp = 1;  % component (Zxx=1, Zxy=2, Zyx=3, Zyy=4).
sidx = 1;  % segment number

figure();
    qopts = copts;
    qqplot_(TFs{3},qopts,comp,fidx);
figure();
    qopts = copts;
    qqplot_(TFs{1},qopts,comp,fidx);



fid = fopen(fullfile(copts.printOptions.printDir,'figures.html'),'w');
fprintf(fid,'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">\n');
fprintf(fid,'<html>\n');
fprintf(fid,'<head>\n');
fprintf(fid,'  <meta http-equiv="Content-type" content="text/html;charset=UTF-8">\n');
fprintf(fid,'  <title>%s</title>\n','title');
fprintf(fid,'</head>\n');

dlist = dir(copts.printOptions.printDir);
[~,idx] = sort([dlist.datenum]);
dlist = dlist(idx);
for i = 1:length(dlist)
    if dlist(i).isdir == 1 || ~endsWith(dlist(i).name,'.png')
        continue;
    end
    fprintf(fid,'<img src="%s" width="500px">\n',dlist(i).name);
end
fclose(fid);
