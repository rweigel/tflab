function [ax1,ax2] = snplot(S,popts)
%SNPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

% Default options
opts = struct();
    opts.title = '';
    opts.print = 0;
    opts.printname = 'snplot';
    opts.printdir = '';
    opts.printfmt = {'pdf'};

    opts.period = 1;
    opts.unwrap = 1;
    opts.vs_period = 1;
    opts.period_range = [];

% Use default options if option not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if iscell(S) && length(S) == 1
    S = S{1};
end
% TODO: Above code is copy of code in zplot().

if isstruct(S)
    % Single transfer function
    figprep();
    fe = S.Metrics.PSD.Smoothed.fe;
    if opts.period
        x = S.Options.info.timedelta./fe;
    else
        x = fe/S.Options.info.timedelta;
    end
    ax1 = subplot('Position', PositionTop);
        for j = 1:size(S.Out,2)
            if iscell(S.Options.info.outstr)
                ls{j} = sprintf('%s\n', S.Options.info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n',S.Options.info.outstr,j);
            end
        end
        semilogx(x,S.Metrics.SN.Smoothed,...
                 'marker','.','markersize',10,'linewidth',2);
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        tflab_title(S,opts,'sn');
        ylabel('Signal to Error');
        grid on;box on;hold on;
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
    ax2 = subplot('Position', PositionBottom);
    
        for j = 1:size(S.In,2)
            if iscell(S.Options.info.outstr)
                ls{j} = sprintf('%s\n', S.Options.info.outstr{j});
            else
                ls{j} = sprintf('%s(:,%d)\n', S.Options.info.outstr,j);
            end
        end
        semilogx(x,S.Metrics.Coherence.Smoothed,...
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
        ylabel('Measured to Predicted Coherence');
        grid on;box on;hold on;
        set(gca,'YLim',[0,1]);
        adjust_exponent('x');
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end

	if ~isempty(S.Options.info.timeunit) && opts.vs_period
        axes(ax1); period_lines();
        axes(ax2); period_lines();
    end
        
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s.%s',opts.printname, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end

if iscell(S)
    % Multiple TFs
    segment_aves = 1;     
    % j = columns of SN (components of output)
    for j = 1:size(S{1}.Metrics.SN.Smoothed,2) 
        if j > 1
            figure();
        end
        figprep();
        set(gcf,'DefaultLegendAutoUpdate','off')
        ax1 = subplot('Position', PositionTop);
            max_T = 0;
            for s = 1:length(S)
                if isfield(S{s}.Metrics,'Segment')
                    segment_aves = 1;
                else
                    segment_aves = 0;
                end
                if segment_aves
                    ye = S{s}.Metrics.Segment.SN_95_boot(:,:,j);
                    SN = mean(S{s}.Metrics.Segment.SN(:,j,:),3);
                    fe = S{s}.Metrics.Segment.fe;
                else
                    SN = S{s}.Metrics.SN.Smoothed(:,j);
                    fe = S{s}.Metrics.PSD.Smoothed.fe;
                    %fe = S{s}.Metrics.fe;
                end
                if opts.period
                    x = 1./fe;
                    max_T = max(max_T,max(x(x < Inf)));
                else
                    x = fe;
                end
                h(s) = semilogx(x,SN,...
                         'marker','.','markersize',10,'linewidth',2);
                grid on;box on;hold on;                 
                ls{s} = sprintf('%s %s',...
                    S{s}.Options.info.stationid,...
                    S{s}.Options.description);
                if segment_aves
                    yl = SN-ye(:,1);
                    yu = -SN+ye(:,2);
                    errorbars(x,SN,yl,yu);
                end

            end                
            legend(h,ls,'Location','NorthEast','Orientation','Horizontal');
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
            adjust_ylim('both');
            limx = get(gca,'XLim');
            plot([limx(1),limx(2)],[1,1],'k');
            if 0
                % export_fig can't handle transparent patches
                patch([limx(1), limx(2), limx(2), limx(1)],...
                      [0, 0, 1, 1],[0.1,0.1,0.1],...
                      'FaceAlpha',0.1,...
                        'LineStyle','none');
            end
            yt = get(gca,'YTick');
            if (yt(1) < 0)
                ytl = get(gca,'YTickLabel');
                ytl{1} = '';
            end
        ax2 = subplot('Position', PositionBottom);
            for s = 1:length(S)
                % TODO: Repeated code
                if isfield(S{s}.Metrics,'Segment')
                    segment_aves = 1;
                else
                    segment_aves = 0;
                end
                if segment_aves
                    ye = S{s}.Metrics.Segment.SN_95_boot(:,:,j);
                    SN = mean(S{s}.Metrics.Segment.SN(:,j,:),3);
                    fe = S{s}.Metrics.Segment.fe;
                else
                    Coh = S{s}.Metrics.Coherence.Smoothed(:,j);
                    fe = S{s}.Metrics.PSD.Smoothed.fe;
                    %fe = S{s}.Metrics.fe;
                end
                if opts.period
                    x = 1./fe;
                    max_T = max(max_T,max(x(x < Inf)));
                else
                    x = fe;
                end
                semilogx(x,Coh,...
                         'marker','.','markersize',10,'linewidth',2);
                grid on;box on;hold on;
            end
            legend(ls,'Location','NorthEast','Orientation','Horizontal');
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
                pre = sprintf('%s(:,%d)',S{1}.Options.info.outstr{j},j);
            end            
            ylabel(sprintf('%s Measured to Predicted Coherence',pre));
            set(gca,'YLim',[0,1]);
            adjust_ylim('both');
            adjust_exponent('x');            
            yt = get(gca,'YTick');
            if (yt(1) < 0)
                ytl = get(gca,'YTickLabel');
                ytl{1} = '';
            end
            
        if ~isempty(S{1}.Options.info.timeunit) && opts.vs_period
            axes(ax1); period_lines();
            axes(ax2); period_lines();
        end
        ext = regexprep(S{1}.Options.info.outstr{j},'\$','');
        if opts.print
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s-%s.%s',opts.printname, ext, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname), opts);
            end
        end
    end
end

