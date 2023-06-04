function [ax1,ax2] = snplot(S,popts)
%SNPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

opts = tflabplot_options(S, popts, '', 'snplot');

if ~iscell(S)
    S = {S};
end

S = defaultinfo(S);
% TODO: Check all same. This assumes 1s.
timeunit = S{1}.Options.info.timeunit;

for s = 1:length(S)
    if 1 || startsWith(opts.type,'averaged')
        fe{s} = S{s}.fe;
        y1{s} = S{s}.Metrics.SN;
        y2{s} = S{s}.Metrics.Coherence;
    end
    if opts.vs_period
        x{s} = S{s}.Options.info.timedelta./fe{s};
    else
        x{s} = fe{s}/S{s}.Options.info.timedelta;
    end
end

if length(S) == 1
    % Single transfer function
    figprep();
    x  = x{1};
    y1 = y1{1};
    y2 = y2{1};
    lg = legend_(S{1});

    ax1 = subplot('Position', opts.PositionTop);
        
        plot(x,y1,opts.line{:});
        colororder_(ax1, y1);
        grid on;box on;hold on;
        if size(y1,2) > 1
            legend(lg,opts.legend{:});
        end
        set(gca,'XTickLabel',[]);
        titlestr(S,opts,'sn');
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if max(y1(:)) > 100
            set(gca,'YScale','log');
        end
        ylabel('Signal to Error');
        yline(1,'k');
        adjust_ylim('upper');
        adjust_yticks();
        adjust_exponent();
        setx(opts,0,timeunit);

    ax2 = subplot('Position', opts.PositionBottom);
        semilogx(x,y2,opts.line{:});
        colororder_(ax2, y2);
        grid on;box on;hold on;
        if size(y2,2) > 1
            legend(lg,opts.legend{:});
        end
        ylabel('Meas. to Pred. Coherence');
        set(gca,'YLim',[0,1]);
        yline(1,'k');
        adjust_ylim('upper');
        adjust_exponent('x');
        setx(opts,1,timeunit);
        
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s.%s',opts.printname, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end

if length(S) > 1
    % Multiple TFs
    comp = 1;
    figprep();
    lg = legend_(S,comp);
    ax1 = subplot('Position', opts.PositionTop);
        grid on;box on;hold on;
        for s = 1:length(y1)
            plot(x{s},y1{s}(:,comp),opts.line{:});
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if max(y1{s}(:,comp)) > 100
            set(gca,'YScale','log');
        end
        legend(lg,opts.legend{:});
        ylabel('Signal to Error');
        adjust_ylim('upper');
        adjust_yticks();
        adjust_exponent();
        yline(1,'k');
        setx(opts,0,timeunit);

    ax2 = subplot('Position', opts.PositionBottom);
        grid on;box on;hold on;
        for s = 1:length(y2)
            plot(x{s},y2{s}(:,comp),opts.line{:});
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        legend(lg,opts.legend{:});
        ylabel('Meas. to Pred. Coherence');
        set(gca,'YLim',[0,1]);
        adjust_ylim('upper');
        adjust_exponent('x');            
        setx(opts,1,timeunit);

    if opts.print
        ext = regexprep(S{1}.Options.info.outstr{comp},'\$','');        
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s-%s.%s',opts.printname, ext, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end
end % function

function lg = legend_(S,comp)

    if iscell(S)
        if (nargin == 1)
            comp = 1;
        end
        for s = 1:length(S)
            desc = S{s}.Options.description;
            if ~isempty(desc)
                desc = [' ',desc];
            end
            if iscell(S{s}.Options.info.outstr)
                lg{s} = sprintf('%s%s\n',...
                            S{s}.Options.info.outstr{comp}, desc);
            else
                if size(S{s}.Out,2) == 1
                    lg{s} = sprintf('%s\n',desc);
                else
                    lg{s} = sprintf('%s(:,%d)%s\n',...
                                S{s}.Options.info.outstr,comp,desc);
                end
            end
        end
    else
        for j = 1:size(S.Out,2)
            if iscell(S.Options.info.outstr)
                lg{j} = sprintf('%s%s\n', S.Options.info.outstr{j});
            else
                if size(S.Out,2) == 1
                    lg{j} = sprintf('%s\n',S.Options.info.outstr);
                else
                    lg{j} = sprintf('%s(:,%d)%s\n',S.Options.info.outstr,j);
                end
            end
        end
    end
    
end


