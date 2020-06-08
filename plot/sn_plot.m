function sn_plot(S,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.period = 1;
    opts.unwrap = 1;
    opts.plottype = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re, Im
    opts.period_range = [];
    opts.savedir = '';
    opts.savefmt = struct();
        opts.savefmt.pdf = 0;

% Use default options if option not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

if length(opts.savedir) > 0 && opts.savedir(end) ~= filesep()
    savedir = [opts.savedir,filesep()];
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if iscell(S) && length(S) == 1
    S = S{1};
end

% TODO: Above code is copy of transferfnZ_plot().

% Single transfer function
if isstruct(S) || length(S) == 1
    figure();
    figprep();
    if opts.period
        x = S.Options.info.timedelta./S.Metrics.fe;
    else
        x = S.Metrics.fe/S.Options.info.timedelta;
    end
    subplot('Position', PositionTop);
        for j = 1:size(S.Out,2)
            if iscell(S.Options.info.outstr)
                ls{j} = sprintf('%s\n', S.Options.info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n',S.Options.info.outstr,j);
            end
        end
        semilogx(x,S.Metrics.SN,...
                 'marker','.','markersize',10,'linewidth',2);
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(S.Options.description,'FontWeight','Normal');        
        ylabel('Signal to Error');
        grid on;box on;hold on;
        if isfield(opts,'title') && ~isempty(opts.title)
            title(opts.title,'FontWeight','normal');
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
    subplot('Position', PositionBottom);
    
        for j = 1:size(S.In,2)
            if iscell(S.Options.info.outstr)
                ls{j} = sprintf('%s\n', S.Options.info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n', S.Options.info.outstr,j);
            end
        end
        semilogx(x,S.Metrics.Coherence,...
                 'marker','.','markersize',10,'linewidth',2);
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('$T$ [%s]', S.Options.info.timeunit));
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        else
            xlabel(sprintf('$f$ [1/%s]', S.Options.info.timeunit));
        end
        ylabel('Coherence');
        grid on;box on;hold on;
        set(gca,'YLim',[0,1]);
        adjust_exponent('x');
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
    figsave(opts,'SN');        
else    
    for j = 1:size(S{1}.Metrics.SN,2)
        figure();
        figprep();
        subplot('Position', PositionTop);
            for s = 1:length(S)
                if opts.period
                    x = 1./S{s}.Metrics.fe;
                else
                    x = S{s}.Metrics.fe;
                end                
                semilogx(x,S{s}.Metrics.SN(:,j),...
                         'marker','.','markersize',10,'linewidth',2);
                grid on;box on;hold on;                 
                ls{s} = S{s}.Options.description;
            end                
            legend(ls,'Location','Best');
            set(gca,'XTickLabel',[]);
            if iscell(S{1}.Options.info.outstr)
                pre = S{1}.Options.info.outstr{j};
            else
                pre = sprintf('%s(:,%d)',S{1}.Options.info.outstr,j);
            end
            ylabel(sprintf('%s Signal to Error',pre));
            if opts.period
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            end            
            if isfield(opts,'title') && ~isempty(opts.title)
                title(opts.title,'FontWeight','normal');
            end
        subplot('Position', PositionBottom);
            for s = 1:length(S)    
                if opts.period
                    x = 1./S{s}.Metrics.fe;
                else
                    x = S{s}.Metrics.fe;
                end                
                semilogx(x,S{s}.Metrics.Coherence(:,j),...
                            'marker','.','markersize',10,'linewidth',2);
                grid on;box on;hold on;
                ls{s} = S{s}.Options.description;
            end
            legend(ls,'Location','Best');
            if opts.period
                xlabel(sprintf('$T$ [%s]', S{1}.Options.info.timeunit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                xlabel(sprintf('$f$ [1/%s]', S{1}.Options.info.timeunit));
            end
            if iscell(S{1}.Options.info.outstr)
                pre = S{1}.Options.info.outstr{j};
            else
                pre = sprintf('%s(:,%d)',S{1}.Options.info.outstr,j);
            end            
            ylabel(sprintf('%s Coherence',pre));
            set(gca,'YLim',[0,1]);
        figsave(opts,
            sprintf('SN-%s',strrep(S{1}.Options.info.outstr,'$','')));            
    end
end