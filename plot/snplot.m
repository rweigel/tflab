function [ax1,ax2] = snplot(S,popts)
%SNPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);

popts = tflabplot_options(S, popts, '', 'snplot');

if ~iscell(S)
    S = {S};
end

% TODO: Check all same. This assumes 1s.
timeunit = S{1}.Metadata.timeunit;

for s = 1:length(S)
    fe{s} = S{s}.Metrics.fe;
    y1{s} = S{s}.Metrics.SN;
    y2{s} = S{s}.Metrics.Coherence;
    if popts.vs_period
        x{s} = S{s}.Metadata.timedelta./fe{s};
    else
        x{s} = fe{s}/S{s}.Metadata.timedelta;
    end
end

if length(S) == 1
    % Single transfer function
    figprep();
    x  = x{1};
    y1 = y1{1};
    y2 = y2{1};
    lg = legend_(S{1},popts);

    ax1 = subplot('Position', popts.PositionTop);
        
        plot(x,y1,popts.line{:});
        colororder_(ax1, y1);
        grid on;box on;hold on;
        if size(y1,2) > 1
            legend(lg,popts.legend{:});
        end
        set(gca,'XTickLabel',[]);
        titlestr(S{1},popts,'sn');
        if popts.vs_period
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
        setx(popts,0,timeunit);

    ax2 = subplot('Position', popts.PositionBottom);
        semilogx(x,y2,popts.line{:});
        colororder_(ax2, y2);
        grid on;box on;hold on;
        if size(y2,2) > 1
            legend(lg,popts.legend{:});
        end
        ylabel('Meas. to Pred. Coherence');
        set(gca,'YLim',[0,1]);
        yline(1,'k');
        adjust_ylim('upper');
        adjust_exponent('x');
        setx(popts,1,timeunit);
        
    if popts.print
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s.%s',popts.printname, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname), popts);
        end
    end
end

if length(S) > 1
    % Multiple TFs
    comp = 1;
    figprep();
    lg = legend_(S,popts,comp);
    ax1 = subplot('Position', popts.PositionTop);
        grid on;box on;hold on;
        for s = 1:length(y1)
            plot(x{s},y1{s}(:,comp),popts.line{:});
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if max(y1{s}(:,comp)) > 100
            set(gca,'YScale','log');
        end
        legend(lg,popts.legend{:});
        ylabel('Signal to Error');
        adjust_ylim('upper');
        adjust_yticks();
        adjust_exponent();
        yline(1,'k');
        setx(popts,0,timeunit);

    ax2 = subplot('Position', popts.PositionBottom);
        grid on;box on;hold on;
        for s = 1:length(y2)
            plot(x{s},y2{s}(:,comp),popts.line{:});
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        legend(lg,popts.legend{:});
        ylabel('Meas. to Pred. Coherence');
        set(gca,'YLim',[0,1]);
        adjust_ylim('upper');
        adjust_exponent('x');            
        setx(popts,1,timeunit);

    if popts.print
        ext = regexprep(S{1}.Options.info.outstr{comp},'\$','');        
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s-%s.%s',popts.printname, ext, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname), popts);
        end
    end
end
end % function

function lg = legend_(S,popts,comp)

    if iscell(S)
        if (nargin == 1)
            comp = 1;
        end
        for s = 1:length(S)
            desc = S{s}.Options.description;
            if ~isempty(desc)
                desc = [' ',desc];
            end
            if iscell(popts.outstr)
                lg{s} = sprintf('%s%s\n',popts.outstr{comp}, desc);
            else
                if size(S{s}.Out,2) == 1
                    lg{s} = sprintf('%s\n',desc);
                else
                    lg{s} = sprintf('%s(:,%d)%s\n',popts.outstr,comp,desc);
                end
            end
        end
    else
        for j = 1:size(S.Out,2)
            if iscell(popts.outstr)
                lg{j} = sprintf('%s%s\n', popts.outstr{j});
            else
                if size(S.Out,2) == 1
                    lg{j} = sprintf('%s\n',popts.outstr);
                else
                    lg{j} = sprintf('%s(:,%d)%s\n',popts.outstr,j);
                end
            end
        end
    end
    
end


