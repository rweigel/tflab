function sn_plot(S1,S2,popts)

% Default options
opts = struct();
    opts.type = 1;
    opts.unwrap = 1;
    opts.magnitude = 0;
    opts.period = 1;
    opts.period_range = [];

% Use default options if option not given
if nargin > 2
    fns = fieldnames(popts);
    for i = 1:length(fns)
        if isfield(opts,fns{i})
           opts.(fns{i}) = popts.(fns{i});
        end
    end
end

figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

info = S1.Options.info;

% Single transfer function
if nargin < 2
    subplot('Position', PositionTop);
    
        for j = 1:size(S1.In,2)
            if iscell(info.outstr)
                ls{j} = sprintf('%s\n', info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n',info.outstr,j);
            end
        end
        if opts.period
            x = 1./S1.fe;
        else
            x = S1.fe;
        end
        semilogx(x,S1.Metrics.SN,...
                 'marker','.','markersize',10,'linewidth',2);
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(S1.Options.description,'FontWeight','Normal');        
        ylabel('Signal to Error');
        grid on;box on;hold on;
        if isfield(popts,'title') && ~isempty(popts.title)
            title(popts.title,'FontWeight','normal');
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
    subplot('Position', PositionBottom);
    
        for j = 1:size(S1.In,2)
            if iscell(info.outstr)
                ls{j} = sprintf('%s\n', info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n',info.outstr,j);
            end
        end
        if opts.period
            x = 1./S1.fe;
        else
            x = S1.fe;
        end
        semilogx(x,S1.Metrics.Coherence,...
                 'marker','.','markersize',10,'linewidth',2);
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('$T$ [%s]', S1.Options.info.timeunit));
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        else
            xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        end
        ylabel('Coherence');
        grid on;box on;hold on;
        set(gca,'YLim',[0,1]);
        adjust_exponent('x');
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
end

if nargin > 1
    
    for j = 1:size(S1.Metrics.SN,2)
        figure();
        figprep();
        if opts.period
            x1 = 1./S1.fe;
            x2 = 1./S2.fe;
        else
            x1 = S1.fe;
            x1 = S2.fe;
        end

        subplot('Position', PositionTop);
    
            semilogx(x1,S1.Metrics.SN(:,j),...
                     'marker','.','markersize',10,'linewidth',2);
            grid on;box on;hold on;                 
            semilogx(x2,S2.Metrics.SN(:,j),...
                     'marker','.','markersize',10,'linewidth',2);

            ls = {S1.Options.description,S2.Options.description};
            legend(ls,'Location','NorthEast');
            set(gca,'XTickLabel',[]);
            if iscell(info.instr)
                pre = info.outstr{j};
            else
                pre = sprintf('%s(:,%d)',info.outstr,j);
            end
            if opts.period
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            end            
            ylabel(sprintf('%s Signal to Error',pre));
            if isfield(popts,'title') && ~isempty(popts.title)
                title(popts.title,'FontWeight','normal');
            end
        subplot('Position', PositionBottom);
    
            semilogx(x1,S1.Metrics.Coherence(:,j),...
                     'marker','.','markersize',10,'linewidth',2);
            grid on;box on;hold on;
            semilogx(x2,S2.Metrics.Coherence(:,j),...
                     'marker','.','markersize',10,'linewidth',2);
            ls = {S1.Options.description,S2.Options.description};
            legend(ls,'Location','NorthEast');
            if opts.period
                xlabel(sprintf('$T$ [%s]', S1.Options.info.timeunit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
            end
            
            if iscell(info.instr)
                pre = info.outstr{j};
            else
                pre = sprintf('%s(:,%d)',info.outstr,j);
            end            
            ylabel(sprintf('%s Coherence',pre));
            grid on;box on;hold on;
            set(gca,'YLim',[0,1]);

    end
end