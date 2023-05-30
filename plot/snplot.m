function [ax1,ax2] = snplot(S,popts)
%SNPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

opts = tflabplot_options(S, popts, '', 'snplot');

if iscell(S) && length(S) == 1
    S = S{1};
end

S = defaultinfo(S);
if iscell(S)
    % TODO: Check all same.
    timeunit = S{1}.Options.info.timeunit;
else
    timeunit = S.Options.info.timeunit;
end

if iscell(S) && length(S) == 1
    S = S{1};
end
% TODO: Above code is copy of code in zplot().

if isstruct(S)
    % Single transfer function
    figprep();
    fe = S.Metrics.PSD.Smoothed.fe;
    if opts.vs_period
        x = S.Options.info.timedelta./fe;
    else
        x = fe/S.Options.info.timedelta;
    end

    lg = legend_(S);

    ax1 = subplot('Position', opts.PositionTop);
        SN = S.Metrics.SN.Smoothed;
        plot(x,SN,opts.line{:});
        colororder_(ax1, SN);
        grid on;box on;hold on;
        if size(S.Out,2) > 1
            legend(lg,opts.legend{:});
        end
        set(gca,'XTickLabel',[]);
        tflab_title(S,opts,'sn');
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if max(SN(:)) > 100
            set(gca,'YScale','log');
        end
        ylabel('Signal to Error');
        yline(1,'k');
        adjust_ylim('upper');
        adjust_yticks();
        adjust_exponent();
        setx(opts,0,timeunit);

    ax2 = subplot('Position', opts.PositionBottom);
        semilogx(x,S.Metrics.Coherence.Smoothed,opts.line{:});
        colororder_(ax2, S.Metrics.Coherence.Smoothed);
        grid on;box on;hold on;
        if size(S.Out,2) > 1
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

if iscell(S)
    % Multiple TFs
    segment_aves = 1;
    % j = columns of SN (components of output)
    for j = 1:size(S{1}.Metrics.SN.Smoothed,2)
        lg = legend_(S, j);
        if j > 1
            figure();
        end
        figprep();
        ax1 = subplot('Position', opts.PositionTop);
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
                if opts.vs_period
                    x = 1./fe;
                    max_T = max(max_T,max(x(x < Inf)));
                else
                    x = fe;
                end
                h(s) = plot(x,SN,opts.line{:});
                grid on;box on;hold on;
                if segment_aves
                    yl = SN-ye(:,1);
                    yu = -SN+ye(:,2);
                    errorbars(x,SN,yl,yu);
                end

            end
            if opts.vs_period
                set(gca,'XScale','log');
            end
            if max(SN(:)) > 100
                set(gca,'YScale','log');
            end
            legend(h,lg,opts.legend{:});
            ylabel('Signal to Error');
            adjust_ylim('upper');
            adjust_yticks();
            adjust_exponent();
            yline(1,'k');
            setx(opts,0,timeunit);
            
        ax2 = subplot('Position', opts.PositionBottom);
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
                if opts.vs_period
                    x = 1./fe;
                    max_T = max(max_T,max(x(x < Inf)));
                else
                    x = fe;
                end
                plot(x,Coh,opts.line{:});
                grid on;box on;hold on;
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
            ext = regexprep(S{1}.Options.info.outstr{j},'\$','');        
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s-%s.%s',opts.printname, ext, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname), opts);
            end
        end
    end
end
end

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


