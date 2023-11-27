function [ax1,ax2] = snplot(S,popts,comp)
%SNPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end
if nargin < 3
    comp = 1;
end

show_xcoh = 1;

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);

popts = tflabplot_options(S,popts,'snplot');

if ~iscell(S)
    S = {S};
end

frequnit = S{1}.Metadata.frequnit;

for s = 1:length(S)
    fe{s} = S{s}.Metrics.fe;
    y1{s} = S{s}.Metrics.SN;
    y2{s} = S{s}.Metrics.Coherence;
    y1cl{s} = S{s}.Metrics.SNCL;
    if show_xcoh == 1
        y3{s} = S{s}.Metrics.Xcoherence;
    end
    if popts.vs_period
        x{s} = 1./(fe{s}*S{s}.Metadata.freqsf);
    else
        x{s} = fe{s}*S{s}.Metadata.freqsf;
    end
end

if length(S) == 1
    % Single transfer function
    figprep();
    x  = x{1};
    y1 = y1{1};
    y1cl = y1cl{1};
    y2 = y2{1};
    y3 = y3{1};
    lg = legend_(S{1},popts);

    ax1 = subplot('Position', popts.PositionTop);
        
        plot(x,y1,popts.line{:});
        colororder_(ax1, y1);
        grid on;box on;hold on;
        if size(y1,2) > 1
            legend(lg,popts.legend{:});
        end
        hold on;
        for i = 1:size(y1,2)
            errorbars(x,y1(:,i),y1(:,i)-y1cl(:,2*i-1),y1cl(:,2*i)-y1(:,i));
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
        setx(popts,0,frequnit);

    ax2 = subplot('Position', popts.PositionBottom);
        semilogx(x,y2,popts.line{:});
        colororder_(ax2, y2);
        grid on;box on;hold on;
        if show_xcoh && length(popts.outstr) == 2 && length(popts.instr) == 2
            % TODO: Generalize; only works for 2x2 Z
            set(gca,'ColorOrderIndex',1);
            h = semilogx(x,y3,'^');
            set(h, {'MarkerFaceColor'}, get(h,'Color'));
            lg = {...
                    '$(E_x,E^{\mathrm{pred}}_x)$','$(E_y,E^{\mathrm{pred}}_y)$',...
                    '$(E_x,B_y)$','$(E_y,B_x)$'...
                  };
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
        
    if popts.print
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s.%s',popts.printname, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname));
        end
    end
end

if length(S) > 1
    % Multiple TFs
    figprep();
    lg = legend_(S,popts,comp);
    ax1 = subplot('Position', popts.PositionTop);
        grid on;box on;hold on;
        for s = 1:length(y1)
            plot(x{s},y1{s}(:,comp),popts.line{:});
        end
        for s = 1:length(y1)
            errorbars(x{s},y1{s}(:,comp),y1{s}(:,comp)-y1cl{s}(:,2*comp-1),y1cl{s}(:,2*comp)-y1{s}(:,comp));
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
        yl = get(gca,'YLim');
        set(gca,'YLim',[0,yl(end)]);
        adjust_exponent();
        yline(1,'k');
        setx(popts,0,frequnit);

    ax2 = subplot('Position', popts.PositionBottom);
        grid on;box on;hold on;
        for s = 1:length(y2)
            plot(x{s},y2{s}(:,comp),popts.line{:});
        end
        if show_xcoh == 1
            lg = legend_(S,popts,comp);
            plot(x{s},y3{s}(:,comp),'ko');
            ylabel(sprintf('Coherence with %s',popts.outstr{comp}));
            if length(popts.instr) == 2
                if comp == 1
                    lg{end+1} = popts.instr{2};
                else
                    lg{end+1} = popts.instr{1};
                end
            end
        else
            ylabel('Meas. to Pred. Coherence');
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        legend(lg,popts.legend{:});
        set(gca,'YLim',[0,1]);
        adjust_ylim('upper');
        adjust_exponent('x');            
        setx(popts,1,frequnit);

    if popts.print
        ext = regexprep(S{1}.Metadata.outstr{comp},'\$','');        
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s-%s.%s',popts.printname, ext, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname));
        end
    end
    if comp < size(S{1}.In,2)
        figure();
        comp = comp + 1;
        ax1 = {ax1};
        ax2 = {ax2};
        [ax1{comp},ax2{comp}] = snplot(S,popts,comp);
    end
end


end % function

function lg = legend_(S,popts,comp)

    if (nargin < 2)
        comp = 1;
    end
    if iscell(S)
        for s = 1:length(S)
            desc = S{s}.Options.description;
            if ~isempty(desc)
                desc = [' ',desc];
            end
            if iscell(popts.outstr)
                lg{s} = sprintf('%s%s',popts.outstr{comp}, desc);
            else
                if size(S{s}.Out,2) == 1
                    lg{s} = sprintf('%s',desc);
                else
                    lg{s} = sprintf('%s(:,%d)%s',popts.outstr,comp,desc);
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


