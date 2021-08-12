function fh = spectrum_plot(S,popts)
%SPECTRUM_PLOT

opts = struct();
    opts.type = 'raw';
    opts.filename = 'spectrum';
    opts.title = '';
    opts.period = 1; % 1 to plot spectra vs period instead of frequency.
    opts.period_range = [];    
    opts.savefmt = {}; % One or more extensions allowed by export_fig
                       % e.g. {'pdf'} or {'svg','png'}.

% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

if ischar(opts.savefmt)
    opts.savefmt = {opts.savefmt};
end

info = S.Options.info;

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

figure();
figprep();

if strcmp(opts.type,'raw')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if isempty(opts.title)
        ts = sprintf('Smoothed PSD of Raw Input%s (top) and Output (bottom)',s1);
    else
        ts = opts.title;
    end
    subplot('Position',PositionTop);
        if opts.period
            loglog(info.timedelta./S.Metrics.fe,S.Metrics.PSD.In,...
                   'marker','.','markersize',10,'linewidth',2);
        else
            loglog(S.Metrics.fe/info.timedelta,S.Metrics.PSD.In,...
                   'marker','.','markersize',10,'linewidth',2);                
        end
        grid on;box on;
        for j = 1:size(S.Metrics.PSD.In,2)
            if iscell(info.instr)
                ls{j} = sprintf('%s [%s]\n',...
                        info.instr{j},info.inunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            info.instr,j,info.inunit);
            end
        end
        if ~isempty(opts.period_range)
            keyboard
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
        adjust_exponent('y');
    subplot('Position',PositionBottom);
        if opts.period
            loglog(info.timedelta./S.Metrics.fe,S.Metrics.PSD.Out,...
                   'marker','.','markersize',10,'linewidth',2);                
        else
            loglog(S.Metrics.fe/info.timedelta,S.Metrics.PSD.Out,...
                   'marker','.','markersize',10,'linewidth',2);                
        end
        grid on;box on;
        for j = 1:size(S.Metrics.PSD.In,2)
            if iscell(info.outstr)
                ls{j} = sprintf('%s [%s]\n',...
                            info.outstr{j},info.outunit);
            else
                ls{j} = sprintf('%s(:,1) [%s]\n',...
                           info.outstr,info.outunit);
            end
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('Period [%s]',info.timeunit));
        else
            xlabel(sprintf('Frequency [1/%s]',info.timeunit));
        end
        adjust_exponent('y');
elseif strcmp(opts.type,'windowed')
    if ~isfield(S,'Window')
        logmsg('Data were not windowed. Not plotting windowed timeseries.\n');
        return;
    end
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if isempty(opts.title)
        ts = sprintf('Smoothed PSD of %s-windowed Input%s (top) and Output (bottom)',...
                        S.Options.td.window.functionstr,s1);
    else
        ts = opts.title;
    end
    subplot('Position',PositionTop);
        if opts.period
            loglog(1./S.fe,S.Window.PSD.In,...
                   'marker','.','markersize',10,'linewidth',2);                
        else
            loglog(S.fe,S.Window.PSD.In,...
                   'marker','.','markersize',10,'linewidth',2);                
        end
        grid on;box on;
        for j = 1:size(S.In,2)
            if iscell(info.instr)
                instr = info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                        instr,info.inunit);
            else
                instr = info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            instr,j,info.inunit);
            end
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        adjust_exponent('y');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        if opts.period
            loglog(1./S.fe,S.Window.PSD.Out,...
                   'marker','.','markersize',10,'linewidth',2);                
        else
            loglog(S.fe,S.Window.PSD.Out,...
                   'marker','.','markersize',10,'linewidth',2);                
        end
        grid on;box on;    
        for j = 1:size(S.Out,2)
            if iscell(info.outstr)
                ls{j} = sprintf('%s [%s]\n',...
                        info.outstr{j},info.outunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            info.outstr,j,info.outunit);
            end
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('Period [%s]',info.timeunit));
        else
            xlabel(sprintf('Frequency [1/%s]',info.timeunit));
        end
elseif strcmp(opts.type,'error')
    ts = sprintf('PSD(error) (top) and PSD(output)./PSD(error) (bottom)');
    subplot('Position',PositionTop);
        loglog(S.Metrics.fe,S.Metrics.PSD.Error);
        grid on;box on;
        if iscell(info.outstr)
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        info.outstr{1},info.outstr{1});
        else
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        info.outstr,info.outstr);
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')        
    subplot('Position',PositionBottom);
        loglog(S.fe,S.Metrics.PSD.Out./S.Metrics.PSD.Error);
        grid on;box on;
        if iscell(info.outstr)
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        info.outstr{1},info.outstr{1});
        else
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        info.outstr,info.outstr);
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        xlabel(sprintf('Frequency [1/%s]',info.timeunit));
end

for i = 1:length(opts.savefmt)
    figsave([opts.filename,'.',opts.savefmt{i}]);
end
