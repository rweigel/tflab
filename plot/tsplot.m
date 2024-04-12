function ax = tsplot(S,popts)
%TSPLOT  Plot timeseries in output of TRANSFERFNFD.
%
%  TSPLOT(S), where S is the output of TRANSFERFNFD or a cell array of
%  such outputs.
%
%  TSPLOT(S, opts) creates a plot using options in structure opts. Use
%  tflabplot_options(S, struct(), 'tsplot') to determine defaults.
%
%  opts.type can be 'original', 'error', 'final', or the name of one of 
%  the fields in S.In_ and S.Out_.

if ischar(S)
    S = loadtf(S);
end

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end
if ~isfield(popts,'type')
    popts.type = 'original';
end

if isfield(popts,'time_range')
    if ~iscell(popts.time_range)
        %popts.time_range = num2cell(popts.time_range);
    end
    try
        fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
        popts.mldatenum_range(1) = datenum(popts.time_range{1},fmt);
        popts.mldatenum_range(2) = datenum(popts.time_range{2},fmt);
    catch
        msg = 'Could not parse one or more elements in popts.time_range.';
        warning('%s Format must be %s',msg,fmt);
    end
end

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);

% Apply default plot options for fields in popts not specified.
popts = tflabplot_options(S, popts, 'tsplot');

tparts = split(popts.type,'-');

final = 0;
if length(tparts) > 1 && strcmp(tparts{2},'final')
    final = 1;
    if ~isfield_(S,'Out_.Final')
        error('Out_.Final (Out data after preprocess) not found .');
    end
    logmsg('Plotting final data.\n');
end

if iscell(S) && length(S) == 1
    S = S{1};
end

figprep();

if iscell(S)
    if ~strcmp(tparts{1},'error')
        error('opts.type for cell array input must be ''error''');
    end
    info = S{1}.Metadata;

    t = index2mldn_(S{1}, size(S{1}.In,1));
    trange = [t(1),t(end)];

    timeunits = {};
    timestarts = {};
    for s = 1:length(S)
        timeunits{s}  = S{s}.Metadata.timeunit;
        timestarts{s} = S{s}.Metadata.timestart;
    end
    % TODO: Allow different timeunits and timestarts.
    if length(unique(timeunits)) > 1
        timeunits
        error('Time units must all be the same');
    end
    if length(unique(timestarts)) > 1
        error('Time starts must all be the same');
    end
else
    info = S.Metadata;
    if strcmp(tparts{1},'error')
        logmsg('Plotting output and predicted for single transfer function.\n');
        if final == 0
            y1{1} = S.Out;
        else
            y1{1} = S.Out_.Final;
        end
        if ~isfield_(S,'Out_.Predicted')
            logmsg('Out_.Predicted and Out_.Error not found. Computing.\n');
            S = tflab_preprocess(S);
            S = tflab_metrics(S);
        end
        y1{2} = S.Out_.Predicted;

        t1 = index2mldn_(S, size(y1{1},1));

        if final == 0
            y2 = S.Out_.Error;
        else
            y2 = S.Out_.ErrorFinal;
        end
        t2 = index2mldn_(S, size(y2,1));

        trange = [t1(1),t1(end)];
        lg1 = {'Measured', 'Predicted'};
    else
        logmsg('Plotting input/output for single transfer function.\n')
        if strcmp(popts.type,'original')
            y1 = S.In;
            y2 = S.Out;
        elseif strcmp(popts.type,'final')
            if ~isfield(S,'In_')
                S = tflab_tdpreprocess(S);
            end
            if ~isfield(S,'In_') || isfield(S.In,'Final')
                msg1 = 'Final data is the same as original b/c no filtering';
                msg2 = ' applied. Request plot for type=original instead.';
                error('%s%s',msg1,msg2);
            end
            y1 = S.In_.Final;
            y2 = S.Out_.Final;
        else
            typeuc = [upper(popts.type(1)),popts.type(2:end)];
            if ~isfield(S,'In_')
                S = tflab_tdpreprocess(S);
            end
            if ~isfield(S,'In_') || ~isfield(S.In_,typeuc)
                msg = 'Invalid plot type request of %s: data were not %s.';
                error(msg, popts.type,popts.type);
            end
            y1 = S.In_.(typeuc);
            y2 = S.Out_.(typeuc);
        end
        t1 = index2mldn_(S, size(y1,1));
        t2 = index2mldn_(S, size(y2,1));
        trange = [t1(1),t1(end)];
        [yl1, yl2] = ylabel_(S,popts);
        [lg1, lg2] = legend_(popts);
    end
end

if ~iscell(S) && ~strcmp(tparts{1},'error')

    ax(1) = subplot('Position',popts.PositionTop);
        plot(t1,y1);
        if size(y1,2) == 1
            % If single line, make black
            colororder(ax(1), {'k'})
        end
        grid on;grid minor;box on;
        titlestr(S,popts,'ts');
        ylabel(yl1);
        if ~isempty(lg1)
            [~, lo] = legend(lg1,popts.legend{:});
            adjust_legend_lines(lo);
        end
        adjust_ylim('upper');
        adjust_exponent('y');
        setx_(0,info,trange);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot(t1,y2);
        if size(y2,2) == 1
            % If single line, make black
            colororder(ax(1), {'k'})
        end
        grid on;grid minor;box on;
        ylabel(yl2);
        if ~isempty(lg2)
            [~, lo] = legend(lg2,popts.legend{:});
            adjust_legend_lines(lo);
        end
        if isfield(S,'OutNoise') && ~strcmp(popts.type,'windowed')
            %hold on;
            %plot(t,S.OutNoise);
        end
        adjust_ylim('upper');
        adjust_exponent('y');
        setx_(1,info,trange);

    figsave_(popts);
end


if ~iscell(S) && strcmp(tparts{1},'error')

    for j = 1:size(S.Out,2)
        if j > 1
            figure();
            figprep();
        end
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        if iscell(popts.outstr)
            outstrerr = popts.outstr{j};
            outstr = sprintf('%s%s',popts.outstr{j},outunit);
        else
            outstr = popts.outstr;
            outstrerr = popts.outstr;
            if size(S.Out,2) > 1
                outstr = sprintf('%s(:,%d)%s',outstr,j,outunit);
                outstrerr = sprintf('%s(:,%d)',outstr,j);
            end
        end
        if final == 0
            Metrics = S.Metrics;
        else
            Metrics = S.MetricsFinal;
        end
        metrics = sprintf('PE/CC/MSE = %.3f/%.3f/%.3f',...
                        Metrics.PE(j),...
                        Metrics.CC(j),...
                        Metrics.MSE(j));

        lg2 = sprintf('%s Error; %s',outstrerr,metrics);

        ax(1,j) = subplot('Position',popts.PositionTop);
            plot(t1,y1{1}(:,j));
            grid on;grid minor;box on;hold on;
            plot(t1,y1{2}(:,j));
            [~, lo] = legend(lg1(:),popts.legend{:});
            ylabel(outstr);
            adjust_legend_lines(lo);
            adjust_ylim();
            adjust_exponent('y');
            setx_(0,info,trange);
            titlestr(S,popts,'ts');
        ax(2,j) = subplot('Position',popts.PositionBottom);
            plot(t2,y2(:,j));
            grid on;grid minor;box on;
            ylabel(outunit);
            %adjust_ylim();
            adjust_exponent('y')
            [~, lo] = legend(lg2,popts.legend{:});
            adjust_legend_lines(lo);
            setx_(1,info,trange);

        figsave_(popts,popts.instr{j});
    end

end

% Compare
if iscell(S)
    if isfield(popts,'mldatenum_range')
        trange = popts.mldatenum_range;
        tidx = find(t >= trange(1) & t <= trange(2));
    else
        tidx = 1:length(t);
    end
    ax(1) = subplot('Position',popts.PositionTop);
        plot(t(tidx),S{1}.Out(tidx,1));
        % Force first line to be black so color order
        % of compared data matches later plots.
        colororder(ax(1), [0,0,0;colororder()]);
        grid on;grid minor;box on;hold on;
        titlestr(S,popts,'ts');
        ylabel(sprintf('%s [%s]',...
                    popts.outstr{1},S{1}.Metadata.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Metadata.stationid);

        for j = 1:length(S)
            plot(t(tidx),S{j}.Out_.Predicted(tidx,1));
            lg{j} = sprintf('Predicted %s\n',...
                        S{j}.Options.description);
        end
        if ~isempty(popts.title)
            title(popts.title);
        end
        [~, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx_(0,info,trange);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot(t(tidx),S{1}.Out(tidx,2));
        % Force first line to be black so color order
        % of compared data matches later plots.
        colororder(ax(2), [0,0,0;colororder()]);
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    popts.outstr{2},S{1}.Metadata.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Metadata.stationid);

        for j = 1:length(S)
            plot(t(tidx),S{j}.Out_.Predicted(tidx,2));
            lg{j} = sprintf('Predicted %s\n',...
                        S{j}.Options.description);
        end

        [~, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx_(1,info,trange);

    if popts.print == 1
        figsave_(popts);
    end
end

end % function

function setx_(last,info,tl)

    if ~isempty(info.timestart)

        % Set tick positions
        set(gca,'XLim',tl);
        adjust_xlim();

        if ischar(info.timestart)
            datetick('x','keeplimits');
            xt = get(gca,'XTick');
            io = 1;
            if xt(1) < tl(1) % First label is not visible.
                io = 2;
            end
            xtl = cellstr(get(gca,'XTickLabel'));
            if xt(2)-xt(1) < 1
                labl = datestr(tl(1),'yyyy/mm/dd');
            else
                labl = datestr(tl(1),'yyyy');
            end
            xtl{io} = sprintf('$$\\begin{array}{c}\\mbox{%s} \\\\ %s\\end{array}$$',...
                                xtl{io},labl);
            set(gca,'XTickLabel',xtl, 'TickLabelInterpreter', 'latex');
        end
    end

    if last == 0
        % Hide tick labels
        set(gca,'XTickLabel',[]);
    end

    if last == 1 && ~ischar(info.timestart)
        % Use default label
        if isempty(info.timeunit)
            xlabel('t');
        else
            xlabel(sprintf('%s since start', info.timeunit));
        end
        adjust_exponent('x');
    end

end

function [lg1, lg2] = legend_(popts)

    lg1 = '';
    lg2 = '';
    if length(popts.instr) > 1
        lg1 = popts.instr;
    end
    if length(popts.outstr) > 1
        lg2 = popts.outstr;
    end
end

function [yl1, yl2] = ylabel_(S,popts)
    meta = S.Metadata;
    if size(S.In,2) > 1
        yl1 = unitstr(meta.inunit);
    else
        yl1 = [popts.instr,' ',unitstr(meta.inunit)];
    end
    if size(S.Out,2) > 1
        yl2 = unitstr(meta.outunit);
    else
        yl2 = [popts.outstr,' ',unitstr(meta.outunit)];
    end
end

function t = index2mldn_(S, nt)

    meta = S.Metadata;
    if ischar(meta.timestart)
        try
            fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
            to = datenum(meta.timestart,fmt);
        catch
            warning(['Could not parse Metadata.timestart. Format must be ',fmt]);
            meta.timestart = '';
        end
        dt = meta.timedelta;
        % Determine ppd (points per day)
        if strcmp(meta.timeunit,'ms') || startsWith(meta.timeunit,'millis')
            ppd = 86400000/dt;
        elseif  startsWith(meta.timeunit,'s')
            ppd = 86400/dt;
        elseif strcmp(meta.timeunit,'m') || startsWith(meta.timeunit,'min')
            ppd = 1440/dt;
        else
            error(['Metadata.timeunit = %s not recognized. ',...
                   'Must be "millseconds", "seconds", or "minutes".'],...
                    meta.timeunit);
        end
        t = to + (0:nt-1)'/ppd;
    else
        t = (meta.timestart:meta.timedelta:nt)';
    end
end

