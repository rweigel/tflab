function timeseries_plot(S,popts)
%TIMESERIES_PLOT - Plot timeseries in output of TRANSFERFNFD.
%
%  TIMESERIES(S, opts), where S is the output of TRANSFERFNFD and pt is the
%  plot type - one of 'raw', 'windowed', or 'error'.

% Default options
opts = struct();
    opts.type = 'raw';
    opts.title = '';
    opts.filename = 'timeseries';
    opts.savefmt = {}; % One or more extensions allowed by export_fig
                       % e.g. {'pdf'} or {'svg','png'}.

% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];
figprep();

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S) && ~isempty(opts.type)
    %warning('opts.type ignored when input is cell array');
    opts.type = '';
end


if iscell(S)
    info = S{1}.Options.info;
    t = S{1}.Time;
    if length(S) > 1
        for j = 2:length(S)
            if ~all(size(S{1}.Time) == size(S{j}.Time))
                error('Time arrays must be identical\n');
            end
            d = S{1}.Time == S{j}.Time;
            if ~all(d(:))
                error('Time arrays must be identical\n');
            end
        end
    end
else
    info = S.Options.info;
    t = S.Time;
end

if ~isempty(info.timestart)
    try
        fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
        to = datenum(info.timestart,fmt);
    catch
        warning(['Could not parse Options.info.timestart. Format must be ',fmt]);
        info.timestart = '';
    end
end

if ~isempty(info.timestart)
    dt = info.timedelta;
    t = [0:length(t)-1]';
    if strcmp(info.timeunit,'ms')
        ppd = 86400000/info.timedelta;
    elseif strcmp(info.timeunit,'s')
        ppd = 86400/info.timedelta;
    elseif strcmp(info.timeunit,'m')
        ppd = 1440/info.timedelta;
    else
        warning('Options.td.timeunit = %s not recognized. Must be ms, s, or m\n',info.timeunit);
        info.timestart = '';
    end
    t = to + t/ppd;
end

if ~iscell(S) && any(strcmp(opts.type,{'raw','windowed','prewhitened'}))

    ts = '';
    if isempty(opts.title)
        site_str = '';
        if isfield(S.Options.info,'stationid')
            sta = S.Options.info.stationid;
            if isempty(sta) == 0
                site_str = sprintf('Site: %s; ',sta);
            end
        end
        if strcmp(opts.type,'raw')
            ts = '';
        end
    else
        ts = opts.title;
    end

    if strcmp(opts.type,'raw')
        In = S.In;
        Out = S.Out;
    end
    
    if strcmp(opts.type,'windowed')
        if ~isfield(S,'Window')
            logmsg('Data were not windowed. Not plotting windowed timeseries.\n');
            return
        end
        if isempty(ts)
            ts = sprintf('%s%s-windowed',...
                         site_str,S.Options.td.window.functionstr);
        end        
        In = S.Window.In;
        Out = S.Window.Out;
    end
    if strcmp(opts.type,'prewhitened')
        if ~isfield(S,'Prewhiten')
            logmsg('Data were not prewhitened. Not plotting prewhitened timeseries.\n');
            return
        end
        if isempty(ts)
            ts = sprintf('%s%s prewhitened',...
                         site_str,S.Options.td.prewhiten.functionstr);
        end        
        In = S.Prewhiten.In;
        Out = S.Prewhiten.Out;
    end
    
    subplot('Position',PositionTop);
    
        plot(t,In);
        grid on;box on;
        ls = legendlabels('In');
        if isfield(S,'InNoise') && ~strcmp(opts.type,'windowed')
            hold on;
            plot(t,S.InNoise);
        end
        title(ts,'FontWeight','Normal');
        [~,h] = legend(ls,'Location','NorthEast','Orientation','Horizontal');
        h = findobj(h,'type','line');
        set(h,'linewidth',2);
        adjust_exponent('y');
        adjust_ylim();
        setx(0,info,[t(1),t(end)]);        
    subplot('Position',PositionBottom);
        plot(t,Out);
        grid on;box on;    
        ls = legendlabels('Out');
        if isfield(S,'OutNoise') && ~strcmp(opts.type,'windowed')
            hold on;
            plot(t,S.OutNoise);
        end
        [~,h] = legend(ls,'Location','NorthEast','Orientation','Horizontal');
        h = findobj(h,'type','line');
        set(h,'linewidth',2);
        adjust_exponent('y');
        adjust_ylim();
        setx(1,info,[t(1),t(end)]);  
    
    for i = 1:length(opts.savefmt)
        figsave([opts.filename,'.',opts.savefmt{i}]);
    end

end

if ~iscell(S) && strcmp(opts.type,'error')

    ts = '';
    if ~isempty(opts.title)
        ts = opts.title;
    end
    
    for j = 1:size(S.Out,2)
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        if iscell(info.outstr)
            outstr = sprintf('%s(:,%d) ',info.outstr{j},j);
        else
            outstr = info.outstr;
        end
        metrics = sprintf('PE/CC/MSE = %.3f/%.3f/%.3f',...
                        S.Metrics.PE(j),...
                        S.Metrics.CC(j),...
                        S.Metrics.MSE(j));
        desc = S.Options.description;
        ls1{j} = sprintf('%s Predicted%s',outstr,outunit);
        ls2{j} = sprintf('%s Error%s; %s',outstr,outunit,metrics);
    end

    subplot('Position',PositionTop);
        plot(t,S.Metrics.Predicted);
        grid on;box on;
        legend(ls1,'Location','NorthEast','Orientation','Vertical');
        adjust_ylim();
        adjust_exponent('y');
        setx(0,info,[t(1),t(end)]);
        title(ts,'FontWeight','Normal');
    subplot('Position',PositionBottom);
        plot(t,S.Metrics.Predicted-S.Out);
        grid on;box on;
        legend(ls2,'Location','NorthEast','Orientation','Vertical');
        adjust_ylim();
        adjust_exponent('y');
        setx(1,info,[t(1),t(end)]);

    filename = [opts.filename,'-A'];
    for i = 1:length(opts.savefmt)
        figsave([filename,'.',opts.savefmt{i}]);
    end

    if size(S.Out,2) == 2
        figure();
        figprep();
        subplot('Position',PositionTop);
            plot(t,S.Out(:,1));
            grid on;box on;hold on;
            plot(t,S.Metrics.Predicted(:,1));        
            legend({ls1{1},ls2{1}},...
                'Location','NorthEast','Orientation','Vertical');
            adjust_ylim();
            adjust_exponent('y');
            setx(0,info,[t(1),t(end)]);
            title(ts);
        subplot('Position',PositionBottom);
            plot(t,S.Out(:,2));
            grid on;box on;hold on;
            plot(t,S.Metrics.Predicted(:,2));
            legend({ls1{2},ls2{2}},...
                    'Location','NorthEast','Orientation','Vertical');
            adjust_ylim();
            adjust_exponent('y');
            setx(1,info,[t(1),t(end)]);

        filename = [opts.filename,'-B'];
        for i = 1:length(opts.savefmt)
            figsave([filename,'.',opts.savefmt{i}]);
        end
    end
end

c = {'k','r','g','b'}; % TODO: Define more colors

% Compare
if iscell(S)
    subplot('Position',PositionTop);
        plot(t,S{1}.Out(:,1),c{1});
        grid on;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    S{1}.Options.info.outstr{1},S{1}.Options.info.outunit));
        ls0 = sprintf('Observed at %s\n', S{1}.Options.info.stationid);

        for j = 1:length(S)
            plot(t,S{j}.Metrics.Predicted(:,1),c{j+1});
            ls{j} = sprintf('Predicted %s\n',...
                        S{j}.Options.description);
        end
        if ~isempty(opts.title)
            title(opts.title,'FontWeight','Normal');
        end
        [lh, lo] = legend({ls0,ls{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(0,info,[t(1),t(end)]);
        %title(ts);
    subplot('Position',PositionBottom);
        plot(t,S{1}.Out(:,2),c{1});
        grid on;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    S{1}.Options.info.outstr{2},S{1}.Options.info.outunit));
        ls0 = sprintf('Observed at %s\n', S{1}.Options.info.stationid);

        for j = 1:length(S)
            for j = 1:length(S)
                plot(t,S{j}.Metrics.Predicted(:,2),c{j+1});
                ls{j} = sprintf('Predicted %s\n',...
                            S{j}.Options.description);
            end
        end
        [lh, lo] = legend({ls0,ls{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(1,info,[t(1),t(end)]);
        
    for i = 1:length(opts.savefmt)
        figsave([opts.filename,'.',opts.savefmt{i}]);
    end
end


function ls = legendlabels(comp)

    unit =  '';
    if strcmp(comp,'In')
        compstrs = info.instr;
        if ~isempty(info.inunit)
            unit = info.inunit;
        end
    else
        compstrs = info.outstr;
        if ~isempty(info.outunit)
            unit = info.outunit;
        end
    end
    unitstr = '';
    if ~isempty(unit)
        unitstr = sprintf(' [%s]',unit);
    end
    for j = 1:size(S.(comp),2)
        if ~iscell(compstrs) && size(S.(comp),2) > 1
            ls{j} = sprintf('%s(:,%d)%s\n',compstrs,j,unitstr);
        else
            if iscell(compstrs)
                ls{j} = sprintf('%s%s\n', compstrs{j}, unitstr);
            else
                ls{j} = sprintf('%s%s\n', compstrs, unitstr);
            end
        end
    end
    
    if isfield(S,[comp,'Noise'])
        jl = j;
        for j = 1:size(S.([comp,'Noise']),2)
            if iscell(compstrs)        
                ls{j+jl} = sprintf('%s Noise%s\n',compstrs{j}, unistr);
            else
                ls{j+jl} = sprintf('%s(:,%d) Noise%s\n',compstrs,j,unitstr);    
            end
        end
    end
end

function setx(last,info,tl)

    if ~isempty(info.timestart)
        % Set tick positions
        set(gca,'XLim',tl);
        adjust_xlim();
        datetick('x','keeplimits');
        xt = get(gca,'XTick');
        io = 1;
        if xt(1) < tl(1) % First label is not visible.
            io = 2; 
        end
        xtl = cellstr(get(gca,'XTickLabel'));
        xtl{io} = sprintf('%s/%s',xtl{io},datestr(tl(1),'yyyy'));
        set(gca,'XTickLabel',xtl);
    end
    if ~last
        % Hide tick labels
        set(gca,'XTickLabel',[]);
    end
    if last && isempty(info.timestart)
        % Use default label
        if ~isempty(info.timeunit)
            xlabel(sprintf('%s since start', info.timeunit));
        else
            xlabel('$t$');
        end
        adjust_exponent();
    end
        
end

end % function

