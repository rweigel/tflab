function [B,E,t,unitsB,unitsE,infile,outfile] = EarthScope_clean(id,plot_,print_)

% List of sites for which cleaning conditions have been created.
prepared = {'VAQ58','ORF03','ORG03'};

if all(strcmp(id,prepared) == 0)
    error('Data for %s is not available',id);
end

if nargin < 2
    plot_ = 0;
end
if nargin < 3
    print_ = 0;
end

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();
addpath(fullfile(scriptdir(),'readers'));

infile  = fullfile(scriptdir(),'data','EarthScope',id,[id,'_raw.mat']);
outfile = strrep(infile,'_raw','_clean');

if exist(outfile,'file') && plot_ == 0 && print_ == 0
    logmsg('Loading: %s\n',relpath(outfile));
    load(outfile);
    logmsg('Loaded:  %s\n',relpath(outfile));
    return;
end

logmsg('Cleaning data for %s\n',id);

data = EarthScope_read(id);

if all(strcmp(id,prepared) == 0)
    error('Data for %s is not available',id);
end

% Create E and B matrices
for c = 1:length(data)
    [datacat{c},timecat{c}] = catSegments(data{c}.segments);
end

B = [datacat{1}, datacat{2}, datacat{3}];
E = [datacat{4}, datacat{5}];
t = timecat{1};

% Number of possible time values.
Ni = 1 + round( (t(end) - t(1))*86400);
logmsg('Number of seconds from t(1) through t(end):      %d\n',Ni);
logmsg('Number of measurements from t(1) through t(end): %d\n',length(t));
% Time value index of measurement
tidx = 1 + round( (t - t(1))*86400 );
t = t(1) + (0:Ni-1)'/86400;
logmsg('Largest gap in seconds:                          %d\n',round(max(diff(tidx))));

% Fill gaps with NaNs.
Ei = nan(Ni,2);
Ei(tidx,:) = E;
E = Ei;

Bi = nan(Ni,3);
Bi(tidx,:) = B;
B = Bi;

if strcmp(id,'VAQ58')
    Ir = (1123950:size(B,1))'; % Values to remove
    Ik = 1:Ir(1)-1;            % Values to keep
    despike_E = {300,[5,5]};   % Despike parameters
    despike_B = {100,[5,5]};   % Despike parameters
end

if strcmp(id,'ORF03')
    Ir = [];
    Ik = [1:size(B,1)];
    %Ir = (1123950:size(B,1))'; % Values to remove
    %Ik = 1:Ir(1)-1;            % Values to keep
    despike_E = {300,[5,5]};   % Despike parameters
    despike_B = {100,[5,5]};   % Despike parameters

    if 0
        delta1 = nanmean(D(850000:end,3))-nanmean(D(11000:848000,3));
        D(848000:end,3) = D(848000:end,3)-delta1;
        D(1:10400,3) = NaN;
        D(852455,5) = NaN;
        D(852480,5) = NaN;
        D(852486,5) = NaN;
        D(861374:861380,5) = NaN;
      end
end

if strcmp(id,'ORG03')
    Ir = [];
    Ik = [1:size(B,1)];
    %Ir = (1123950:size(B,1))'; % Values to remove
    %Ik = 1:Ir(1)-1;            % Values to keep
    despike_E = {300,[5,5]};   % Despike parameters
    despike_B = {100,[5,5]};   % Despike parameters
end

unitsB = 'Counts';
if isfield(data{1}.segments(1),'dataScaledUnits');
    % Assumes units same for all segments & B channels
    unitsB = data{1}.segments(1).dataScaledUnits;
    if strcmp(lower(unitsB),'t')
        B = B/1e-9;
        unitsB = 'nT';
    end
end

unitsE = 'Counts';
if isfield(data{1}.segments(1),'dataScaledUnits');
    % Assumes units same for all segments & E channels
    unitsE = data{4}.segments(1).dataScaledUnits;
    if strcmp(lower(unitsE),'v/m')
        E = E/1e-6;
        unitsE = 'mV/km';
    end
    if strcmp(lower(unitsE),'v')
        warning('Check E units; V given, V/m assumed');
        E = E/1e-6;
        unitsE = 'mV/km';
    end
end

logmsg('Cleaning E');
E1 = E;
E1(Ir,:) = NaN;
E2 = removemean(E1);
E3 = despike(E2, despike_E{:});
E4 = naninterp1(E3);

logmsg('Cleaning B');
B1 = B;
B1(Ir,:) = NaN;
B2 = removemean(B1);
B3 = despike(B2, despike_B{:});
B4 = naninterp1(B3);

to = t;
[B, E, t] = trimnans(B4, E4, t);

infile = relpath(infile);
outfileo = outfile;
outfile = relpath(outfile);
logmsg('Saving: %s\n',relpath(outfile));
save(outfileo,'B','E','t','unitsB','unitsE','infile','outfile');
logmsg('Saved:  %s\n',relpath(outfile));

if plot_ || print_
    figure(1);figprep();clf;

    subplot(4,1,1)
        plot(to,E1);
        datetick();
        grid minor;
        title('Raw $\mathbf{E}$');
        grid on;
        ylabel(unitsE);
        axis tight;
        xlim = get(gca,'XLim');
        legend('$E_x$','$E_y$');
        zoominfo('E raw');

    subplot(4,1,2)
        plot(to,E2)
        datetick();
        grid minor;
        title('Mean subtracted');
        grid on;
        ylabel(unitsE);
        axis tight;
        set(gca,'XLim',xlim);
        legend('$E_x$','$E_y$');
        zoominfo('E after mean subtraction');

    subplot(4,1,3)
        plot(to,E3);
        datetick();
        grid minor;
        title('Despiked and chunk(s) removed');
        grid on;
        ylabel(unitsE);
        axis tight;
        set(gca,'XLim',xlim);
        legend('$E_x$','$E_y$');
        zoominfo('E after despiking and chunk removal');

    subplot(4,1,4)
        plot(to,E4)
        datetick();
        grid minor;
        title('Interpolated over NaNs');
        grid on;
        ylabel(unitsE);
        axis tight;
        set(gca,'XLim',xlim);
        legend('$E_x$','$E_y$');
        zoominfo('E after interpolating over NaNs');
        xlabel(['Month/Day of ',datestr(t(1),'yyyy')])
end

if plot_
    figure(2);figprep();clf;

    subplot(4,1,1)
        plot(to,B1(:,1:2))
        datetick();
        grid minor;
        title('Raw $\mathbf{B}$');
        grid on;
        ylabel(unitsB);
        axis tight;
        legend('$B_x$','$B_y$');
        zoominfo('B');
        xlim = get(gca,'XLim');

    subplot(4,1,2)
        plot(to,B2(:,1:2))
        datetick();
        grid minor;
        title('Mean subtracted');
        grid on;
        ylabel(unitsB);
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');
        zoominfo('B after mean subtraction');

    subplot(4,1,3)
        plot(to,B3(:,1:2));
        datetick();
        grid minor;
        title('Desiked and chunk(s) removed');
        ylabel(unitsB);
        grid on;
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');
        zoominfo('B after despiking and chunk removal');

    subplot(4,1,4)
        plot(to,B4(:,1:2))
        datetick();
        grid minor;
        title('Interpolated');
        grid on;
        ylabel(unitsB);
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');
        xlabel(['Month/Day of ',datestr(t(1),'yyyy')])
        zoominfo('B after interpolation');
end

if print_
    figure(2);
    fname = sprintf('%s_B.png',outfile(1:end-4));
    logmsg(sprintf('Writing: %s\n',relpath(fname)));
	print(fname,'-dpng','-r300');
    logmsg(sprintf('Wrote:   %s\n',relpath(fname)));

    figure(1);
    fname = sprintf('%s_E.png',outfile(1:end-4));
    logmsg(sprintf('Writing: %s\n',relpath(fname)));
	print(fname,'-dpng','-r300');
    logmsg(sprintf('Wrote:   %s\n',relpath(fname)));
end

function [data, time] = catSegments(segments)

    % Does not check that units the same for all segments.
    data = [];
    time = [];
    for s = 1:length(segments)
        logmsg('Channel %s, segment %d: %s to %s\n',...
            segments(s).channel,...
            s,...
            datestr(datevec(segments(s).startTime),31),...
            datestr(datevec(segments(s).endTime),31));

        % Segment data
        if isfield(segments(s),'dataScaled')
            data = [data; segments(s).dataScaled];
        else
            data = [data; segments(s).data];
        end

        % Segment time in fraction of day since start of segment.
        stime = (0:length(segments(s).data)-1)'*(segments(s).sampleRate/86400);

        % Segment time as MATLAB datenum.
        stime = segments(s).startTime + stime;
        time = [time; stime];
    end
end

function zoominfo(varstr)
    zoom off;
    hB = zoom(gca);
    set(hB,'ActionPreCallback',  @msgpre);
    set(hB,'ActionPostCallback', @msgpost);
    function msgpost(obj,evd)
        xlims = round(evd.Axes.XLim);
        minx = max(0,xlims(1));
        fprintf('Showing %s over XLim = [%d,%d]\n', varstr, minx, xlims(2));
    end
    function msgpre(obj,evd), end
end

end

