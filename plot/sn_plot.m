function sn_plot(S,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.period = 1;
    opts.unwrap = 1;
    opts.plottype = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re, Im
    opts.period_range = [];
    opts.savedir = '';
    opts.filename = 'SN';
    opts.savefmt = {};

% Use default options if option not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

if length(opts.savedir) > 0 && opts.savedir(end) ~= filesep()
    opts.savedir = [opts.savedir,filesep()];
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if iscell(S) && length(S) == 1
    S = S{1};
end

% TODO: Above code is copy of code in transferfnZ_plot().

if isstruct(S) % Single transfer function
    ts = opts.title;
    if isempty(ts)
        sta = '';
        if isfield(S.Options.info,'stationid')
            sta = S.Options.info.stationid;
        end
        ts = sprintf('Site: %s',sta);
    end
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
        title(ts,'FontWeight','Normal');        
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
    for i = 1:length(opts.savefmt)
        figsave([opts.filename,'.',opts.savefmt{i}]);
    end
else % Multiple TFs
    segment_aves = 1;
     
    % j = columns of SN (components of output)
    for j = 1:size(S{1}.Metrics.SN,2) 
        figure();
        figprep();
        set(gcf,'DefaultLegendAutoUpdate','off')
        p1 = subplot('Position', PositionTop);
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
                    SN = S{s}.Metrics.SN(:,j);
                    fe = S{s}.Metrics.fe;
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
        p2 = subplot('Position', PositionBottom);
            for s = 1:length(S)    
                if segment_aves
                    Metrics = S{s}.Segment.Metrics;
                else
                    Metrics = S{s}.Metrics;
                end
                if opts.period
                    x = 1./Metrics.fe;
                else
                    x = Metrics.fe;
                end                
                semilogx(x,Metrics.Coherence(:,j),...
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
            ylabel(sprintf('%s Coherence',pre));
            set(gca,'YLim',[0,1]);
            adjust_ylim('both');

        axes(p1); period_lines(max_T);
        axes(p2); period_lines(max_T);            
        ext = regexprep(S{1}.Options.info.outstr{j},'\$','');
        for i = 1:length(opts.savefmt)
            figsave([opts.filename,'-',ext,'.',opts.savefmt{i}]);
        end
    end
end