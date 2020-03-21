function timeseries_plot(S,popts)
%TIMESERIES_PLOT - Plot timeseries in output of TRANSFERFNFD.
%
%  TIMESERIES(S, opts), where S is the output of TRANSFERFNFD and pt is the
%  plot type - one of 'raw', 'windowed', or 'error'.

% Default options
opts = struct();
    opts.type = 'raw';
    opts.title = '';

% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

figure();
figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

info = S.Options.info;
t = S.Time;

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

if size(S.In,2) > 1,
    s1 = 's';
else
    s1 = '';
end

if strcmp(opts.type,'raw') || strcmp(opts.type,'windowed')
    
    if isempty(opts.title)
        ts = sprintf('Raw Input%s (top) and Raw Output (bottom)', s1);
        if strcmp(opts.type,'windowed')
            ts = sprintf('%s-windowed Input%s (top) and Output (bottom)',...
                         S.Options.td.window.functionstr,...
                         s1);
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
        In = S.Window.In;
        Out = S.Window.Out;
    end
    
    subplot('Position',PositionTop);
        plot(t,In);
        grid on;box on;
        for j = 1:size(In,2)
            if iscell(info.instr)        
                ls{j} = sprintf('%s [%s]\n', info.instr{j}, info.inunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',info.instr,j,info.inunit);
            end
        end
        title(ts,'FontWeight','Normal');
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_exponent('y');
        adjust_ylim();
        setx(0,info,[t(1),t(end)]);        
    subplot('Position',PositionBottom);
        plot(t,Out);
        grid on;box on;    
        for j = 1:size(Out,2)
            if iscell(info.outstr)        
                ls{j} = sprintf('%s [%s]\n', info.outstr{j}, info.outunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',info.outstr,j,info.outunit);    
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_exponent('y');
        adjust_ylim();
        setx(1,info,[t(1),t(end)]);
        
elseif strcmp(opts.type,'error')

    if isempty(opts.title)
        ts = sprintf('Method: %s',S.Options.description);
    else
        ts = opts.title;
    end    

    subplot('Position',PositionTop);
        plot(t,S.Out);
        grid on;box on;
        for j = 1:size(S.Out,2)
            if iscell(info.outstr)        
                ls{j} = sprintf('%s [%s]\n', info.outstr{j}, info.outunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',info.outstr,j,info.outunit);    
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Vertical');
        adjust_ylim();
        adjust_exponent('y');
        set(gca,'XTickLabel',[]);
        setx(0,info,[t(1),t(end)]);
        title(ts);
    subplot('Position',PositionBottom);
        plot(t,S.Out-S.Metrics.Predicted);
        grid on;box on;
        for j = 1:size(S.Out,2)        
            metrics = sprintf('PE/CC/MSE = %.3f/%.3f/%.3f',...
                            S.Metrics.PE(j),...
                            S.Metrics.CC(j),...
                            S.Metrics.MSE(j));
            if iscell(info.outstr)
                ls{j} = sprintf('(%s observed) - (%s predicted) [%s]; %s',...
                            info.outstr{j},...
                            info.outstr{j},...
                            info.outunit,...
                            metrics);
            else
                ls{j} = sprintf('(%s(:,1) observed) - (%s(:,1) predicted) [%s]; %s',...
                            info.outstr,...
                            info.outstr,...
                            info.outunit,...
                            metrics);    
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Vertical');
        adjust_ylim();
        adjust_exponent('y');
        setx(1,info,[t(1),t(end)]);
end

drawnow;

figsave(opts,['timeseries-',opts.type]);

function setx(last,info,tl)

    if ~isempty(info.timestart)
        set(gca,'XLim',tl);
        datetick('x','keeplimits');
    end
    if ~last
        set(gca,'XTickLabel',[]);
    else
        xlabel(sprintf('%s since start', info.timeunit));
        adjust_exponent();
    end
        
end

end % function

