function [ax1,ax2] = snplot(S,popts,comps)
%SNPLOT
%
%   SNPLOT(S)
%   SNPLOT(S,opts)
%   SNPLOT(S,opts,component)

msg = 'S must be a tflab struct or cell array of tflab structs';
assert(isstruct(S) || iscell(S), msg);

if isstruct(S)
    S = {S};
end

if nargin < 2
    popts = struct();
end
if nargin < 3
    comps = 1:size(S{1}.Metrics.SN,2);
end
comps = sort(comps);
comp = comps(1);

if length(comps) > 1
    for c = 1:length(comps)
        if c > 1
            figure;
        end
        logmsg('Plotting component %d.\n',comps(c))
        [ax1(c),ax2(c)] = snplot(S,popts,comps(c));
    end
    return
end

show_xcoh = 1;

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);
popts = tflabplot_options(S,popts,'snplot');

frequnit = S{1}.Metadata.frequnit;

onfinal = 0;

if length(S) == 1

    % Single transfer function
    figprep();
    lg = legend_(S{1},popts);

    if (nargin > 2)
        lg = lg{comp};
    end

    [x,y1,y2,y3,y1clu,y1cll] = xyvals_(S,popts,onfinal);
    if nargin > 2
        y1 = y1(:,comp);
        y1cll = y1cll(:,comp);
        y1clu = y1clu(:,comp);
        y2 = y2(:,comp);
    end

    ax1 = subplot('Position', popts.PositionTop);

        plot(x,y1,popts.line{:});
        colororder_(ax1, y1);
        grid on;box on;hold on;
        legend(lg,popts.legend{:});
        hold on;
        for j = 1:size(y1,2)
            errorbar_(x,y1(:,j),y1(:,j)-y1cll(:,j),y1clu(:,j)-y1(:,j));
        end
        set(gca,'XTickLabel',[]);
        title_(S{1},popts,'sn');
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
        setx(popts,0,frequnit);

    ax2 = subplot('Position', popts.PositionBottom);
        semilogx(x,y2,popts.line{:});
        colororder_(ax2, y2);
        grid on;box on;hold on;

        if show_xcoh
            if size(y3,2) == 1
                h = semilogx(x,y3(:,1),'rs');
                set(h, 'MarkerFaceColor', get(h,'Color'));
                instr = replace(popts.instr{1},'$','');
                outstr = replace(popts.outstr{1},'$','');
                a = sprintf('$(%s,%s^{\\mathrm{pred}})$',outstr,outstr);
                b = sprintf('$(%s,%s)$',outstr,instr);
                lg = {a, b};
            end
            if size(y3,2) == 4 && nargin > 2
                % Plotting a single component
                %set(gca,'ColorOrderIndex',1);
                if comp == 1
                    idx = 1:2;
                    lg = {...
                            '$(E_x,E^{\mathrm{pred}}_x)$',...
                            '$(E_x,B_x)$',...
                            '$(E_x,B_y)$'...
                        };
                end
                if comp == 2
                    idx = 3:4;
                    lg = {...
                            '$(E_y,E^{\mathrm{pred}}_y)$',...
                            '$(E_y,B_x)$',...
                            '$(E_y,B_y)$'...
                        };
                end
                h = semilogx(x,y3(:,idx(1)),'rs');
                set(h, 'MarkerFaceColor', get(h,'Color'));
                h = semilogx(x,y3(:,idx(2)),'b^');
                set(h, 'MarkerFaceColor', get(h,'Color'));
            end
            legend(lg,popts.legend{:});
            ylabel('Coherence');
        else
            if size(y2,2) > 1
                legend(lg,popts.legend{:});
            end
            ylabel('Meas. to Pred. Coherence');
        end
        set(gca,'YLim',[0,1]);
        yline(1,'k');
        adjust_ylim('upper');
        adjust_exponent('x');
        setx(popts,1,frequnit);

    figsave_(popts,S{1}.Metadata.outstr{comp})
end

if length(S) > 1
    % Multiple TFs
    figprep();

    lg = legend_(S,popts);
    [x,y1,y2,y3,y1clu,y1cll] = xyvals_(S,popts,onfinal);

    if iscell(popts.outstr)
        outstr = popts.outstr{comp};
    else
        if size(S.Out,2) == 1
            outstr = sprintf('%s\n',popts.outstr);
        else
            outstr = sprintf('%s(:,%d)%s\n',popts.outstr,j);
        end
    end

    ax1 = subplot('Position', popts.PositionTop);
        grid on;box on;hold on;
        c = colororder;
        for s = 1:length(y1)
            marker = markeropts(length(S),s);
            plot(x{s},y1{s}(:,comp),marker{:},'color', c(s,:));
        end
        for s = 1:length(y1)
            plot(x{s}, y2{s}(:,comp)./(1-y2{s}(:,comp)), '-','color', c(s,:));
        end
        title_(S,popts,'sn');
        for s = 1:length(y1)
            errorbar_(x{s},y1{s}(:,comp),y1{s}(:,comp)-y1cll{s}(:,comp),y1clu{s}(:,comp)-y1{s}(:,comp));
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if max(y1{s}(:,comp)) > 100
            set(gca,'YScale','log');
        end
        legend(lg,popts.legend{:});
        ylabel(sprintf('%s Signal to Error',outstr));
        adjust_ylim('upper');
        adjust_yticks();
        yl = get(gca,'YLim');
        set(gca,'YLim',[0,yl(end)]);
        adjust_exponent();
        yline(1,'k');
        setx(popts,0,frequnit);
    ax2 = subplot('Position', popts.PositionBottom);
        grid on;box on;hold on;
        for s = 1:length(y2)
            marker = markeropts(length(S),s);
            plot(x{s},y2{s}(:,comp),marker{:});
        end
        if show_xcoh == 1
            lg = legend_(S,popts);
            for i = 1:length(lg)
                lg{i} = sprintf('%s %s',popts.outstr{comp},lg{i});
            end
            if length(popts.instr) == 2
                if comp == 1
                    plot(x{1},y3{1}(:,2),'ko');
                    lg{end+1} = popts.instr{2};
                else
                    plot(x{1},y3{1}(:,4),'ko');
                    lg{end+1} = popts.instr{1};
                end
            end
            ylabel(sprintf('Squared Coherence with %s',outstr));
        else
            tmp = '%s Meas. to %s Pred. Coherence';
            ylabel(sprintf(tmp,outstr,outstr));
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        legend(lg,popts.legend{:});
        set(gca,'YLim',[0,1]);
        adjust_ylim('upper');
        adjust_exponent('x');
        setx(popts,1,frequnit);

    figsave_(popts,S{1}.Metadata.outstr{comp})

end

end % function

function [x,y1,y2,y3,y1cll,y1clu] = xyvals_(S,popts,onfinal)

    if onfinal == 0
        MetricsName = 'Metrics';
    else
        MetricsName = 'MetricsFinal';
    end

    for s = 1:length(S)
        fe{s} = S{s}.(MetricsName).fe;
        y1{s} = S{s}.(MetricsName).SN;
        y2{s} = sqrt(S{s}.(MetricsName).PredictionCoherence);
        y1clu{s} = S{s}.(MetricsName).SNCLu;
        y1cll{s} = S{s}.(MetricsName).SNCLl;
        y3{s} = sqrt(S{s}.(MetricsName).CrossCoherence);
        if popts.vs_period
            x{s} = 1./(fe{s}*S{s}.Metadata.freqsf);
        else
            x{s} = fe{s}*S{s}.Metadata.freqsf;
        end

        if popts.vs_period && ~isempty(popts.period_range)
            idx = x{s} >= popts.period_range(1) & x{s} <= popts.period_range(2);
            x{s} = x{s}(idx);
            y1{s} = y1{s}(idx,:);
            y2{s} = y2{s}(idx,:);
            y3{s} = y3{s}(idx,:);
            y1clu{s} = y1clu{s}(idx,:);
            y1cll{s} = y1cll{s}(idx,:);
        end
    end

    if length(S) == 1
        x  = x{1};
        y1 = y1{1};
        y1cll = y1cll{1};
        y1clu = y1clu{1};
        y2 = y2{1};
        y3 = y3{1};
    end
end

function lg = legend_(S,popts)

    if iscell(S)
        for s = 1:length(S)
            desc = S{s}.Options.description;
            lg{s} = sprintf('%s',desc);
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
