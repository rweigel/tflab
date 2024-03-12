function [ax1,ax2] = zplot(S,popts,comps)
%ZPLOT
%
%   ZPLOT(S)
%   ZPLOT(S, popts)
%   ZPLOT(S, popts, components)
%
%   opts.type is one of 1, 2, or 3, where
%
%   1 => Z,phase
%   2 => rho,phase (assumes Z in (mV/km)/nT and frequencies in Hz).
%   3 => Real,Imaginary

msg = 'S must be a tflab struct or cell array of tflab structs';
assert(isstruct(S) || iscell(S), msg);

if isstruct(S)
    S = {S};
end

if nargin < 2
    popts = struct();
end
if nargin < 3
    comps = 1:size(S{1}.Z,2);
end
comps = sort(comps);

S = tflab_metadata(S);

% TODO: Check all same units and sfs or allow different.
frequnit = S{1}.Metadata.frequnit;
popts   = tflabplot_options(S, popts, 'zplot');
Zstrs   = popts.zstrs;
Rhostrs = popts.rhostrs;
Phistrs = popts.phistrs;

% Single transfer function
if length(S) == 1
    logmsg('Plotting single transfer function.\n')

    yl1 = unitstr_(S{1}.Metadata);
    yl2 = '[$^\circ$]';
    if popts.type == 1 || popts.type == 2
        template2 = '$%s$';
        strings2 = Phistrs;
        if popts.type == 1
            template1 = '$|%s|$';
            strings1 = Zstrs;
        else
            template1 = '$%s$';
            strings1 = Rhostrs;
        end
    else
        template1 = 'Re$(%s)$';
        strings1 = Zstrs;
        template2 = 'Im$(%s)$';
        strings2 = Zstrs;
        yl2 = unitstr_(S{1}.Metadata);
    end

    for j = 1:length(comps)
        ls1{j} = sprintf(template1,strings1{comps(j)});
        [x1,y1(:,j),dyu1(:,j),dyl1(:,j)] = xyvals_(S,popts,comps(j),'top','parametric');
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            [~,~,dyu1b(:,j),dyl1b(:,j)] = xyvals_(S,popts,comps(j),'top','bootstrap');
        end

        ls2{j} = sprintf(template2,strings2{comps(j)});
        [x2,y2(:,j),dyu2(:,j),dyl2(:,j)] = xyvals_(S,popts,comps(j),'bottom','parametric');
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            [~,~,dyu2b(:,j),dyl2b(:,j)] = xyvals_(S,popts,comps(j),'bottom','bootstrap');
        end
    end

    figprep();
    ax1 = subplot('Position', popts.PositionTop);

        popts.line = markeropts(size(S{1}.Z,1), 1);

        h1 = plot(x1, y1, popts.line{:});
        hold on;grid on;box on;
        titlestr(S{1},popts,'z');
        colororder_(ax1,y1);

        if length(ls1) > 1
            ylabel(yl1);
            legend(ls1,popts.legend{:});
        else
            ylabel(sprintf('%s %s',ls1{1},yl1));
        end

        if popts.vs_period
            set(gca,'XScale','log');
        end

        if popts.type < 3
            set(gca,'YScale','log');
        end

        % Must be called after scale type is set.
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            plot_errorbars_(h1,0.97*x1,y1,dyl1,dyu1,1);
            plot_errorbars_(h1,x1*1.03,y1,dyl1b,dyu1b,1);
        else
            plot_errorbars_(h1,x1,y1,dyl1,dyu1,2);
        end

        adjust_yticks(1e-4);
        adjust_exponent('y');
        setx(popts,0,frequnit);

        if popts.type == 3
            xlims = get(gca(),'XLim');
            line(xlims,[0,0])
        end

    ax2 = subplot('Position', popts.PositionBottom);

        h2 = plot(x2, y2, popts.line{:}, 'markersize', 20);

        hold on;grid on;box on;
        colororder_(ax2,y2);

        if length(ls2) > 1
            ylabel(yl2);
            legend(ls2,popts.legend{:});
        else
            ylabel(sprintf('%s %s',ls2{1},yl2));
        end

        if popts.vs_period
            set(gca,'XScale','log');
        end

        if ~popts.unwrap && popts.type ~= 3
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        end

        % Must be called after scale type is set.
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            plot_errorbars_(h2,0.97*x2,y2,dyl2,dyu2,1);
            plot_errorbars_(h2,x2*1.03,y2,dyl2b,dyu2b,1);
        else
            plot_errorbars_(h2,x2,y2,dyl2,dyu2,2);
        end

        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,1,frequnit);

        xlims = get(gca(),'XLim');
        line(xlims,[0,0])

    if popts.print == 1
        exts = {};
        if length(Zstrs) == 2
            exts = {'x','y'};
        end
        if length(Zstrs) == 4
            exts = {'xx','xy','yx','yy'};
        end
        ext = '';
        for i = 1:length(comps)
            if ~isempty(exts)
                ext = [ext,exts{comps(i)},'-'];
            else
                ext = [ext,num2str(comps(i)),'-'];
            end
        end
        figsave_(popts,ext);
    end

end % if isstruct(S)

% Multiple transfer functions
if length(S) > 1
    logmsg('Plotting multiple transfer functions.\n')

    if length(comps) == 1
        comp = comps;
    else
        for c = 1:length(comps)
            if c > 1
                figure;
            end
            logmsg('Plotting component %d.\n',comps(c))
            [ax1(c),ax2(c)] = zplot(S,popts,comps(c));
        end
        return
    end

    figprep();

    ax1 = subplot('Position', popts.PositionTop);

        [x,y,dyu,dyl] = xyvals_(S,popts,comp,'top','parametric');

        for s = 1:length(S)
            ls{s} = sprintf('%s',S{s}.Options.description);
            popts.line = markeropts(size(S{s}.Z,1), s);
            h(s) = plot(x{s}, y{s}, popts.line{:});
            if s == 1
                grid on;hold on;box on;
            end
        end

        legend(h,ls,popts.legend{:});

        if popts.type < 3
            set(gca,'YScale','log');
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        yunitstr_ = unitstr_(S{1}.Metadata);
        if popts.type == 1
            yl = sprintf('$|%s|$%s',Zstrs{comp},yunitstr_);
        end
        if popts.type == 2
            yl = sprintf('$%s$ [$\\Omega\\cdot$m]',Rhostrs{comp});
        end
        if popts.type == 3
            yl = sprintf('Re$(%s)$ %s',Zstrs{comp},yunitstr_);
        end
        if ~isempty(popts.period_range)
            set(gca,'XLim',popts.period_range);
        end
        ylabel(yl);

        if length(x) == 2
            % Must be called after scale type is set.
            plot_errorbars_(h(1),0.97*x{1},y{1},dyl{1},dyu{1},1);
            plot_errorbars_(h(2),1.03*x{2},y{2},dyl{2},dyu{2},1);
        end

        adjust_yticks(1e-4);
        adjust_exponent('y');

        if popts.type == 3
            adjust_ylim('both');
        else
            adjust_ylim('upper');
        end
        setx(popts,0,frequnit);

    ax2 = subplot('Position', popts.PositionBottom);

        [x,y,dyu,dyl] = xyvals_(S,popts,comp,'bottom','parametric');

        for s = 1:length(S)
            ls{s} = sprintf('%s',S{s}.Options.description);

            popts.line = markeropts(size(S{s}.Z,1), s);
            if popts.vs_period
                h(s) = semilogx(x{s}, y{s}, popts.line{:});
            else
                h(s) = plot(x{s}, y{s}, popts.line{:});
            end

            if s == 1,grid on;box on;hold on;end
        end

        yl = sprintf('$%s$ [$^\\circ]$',Phistrs{comp});
        if popts.type == 3
            yl = sprintf('Im$(%s)$ %s',Zstrs{comp},unitstr_(S{1}.Metadata));
        end
        ylabel(yl);
        %legend(ls,popts.legend{:});

        if popts.type ~= 3 && ~popts.unwrap
            set(gca,'YScale','linear');
            set(gca,'YLim',[-180,180]);
            set(gca,'YTick',-180:45:180);
            adjust_ylim();
        else
            adjust_ylim('upper');
        end

        if length(x) == 2
            % Must be called after scale type is set.
            plot_errorbars_(h(1),0.97*x{1},y{1},dyl{1},dyu{1},1);
            plot_errorbars_(h(2),1.03*x{2},y{2},dyl{2},dyu{2},1);
        end

        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,1,frequnit);

    figsave_(popts,Zstrs{comp});

end % if iscell(S)

end % function zplot()

function plot_errorbars_(h,x,y,dyl,dyu,lw)
    if isempty(dyl),return,end
    colors = get(h,'Color');
    if ~iscell(colors)
        colors = {colors};
    end
    for comp = 1:size(y,2)
        args = {'y','Color',colors{comp},'LineWidth',lw};
        errorbar_(x,y(:,comp),dyl(:,comp),dyu(:,comp),args{:});
    end
end

function [x,y,dyu,dyl] = xyvals_(S,popts,comp,panel,errorbar_method)

    % TODO: Assumes units are the same for all Zs.
    %       Check this before plotting.

    for s = 1:length(S)

        if strcmp(errorbar_method,'parametric')
            ErrorEstimates = S{s}.Metrics.ErrorEstimates.Parametric;
        end
        if strcmp(errorbar_method,'bootstrap')
            ErrorEstimates = S{s}.Metrics.ErrorEstimates.Bootstrap;
        end
        if strcmp(panel,'top')
            if popts.type == 1
                y{s} = abs(S{s}.Z(:,comp));
                dyl{s} = y{s} - ErrorEstimates.ZMAGCL95l(:,comp);
                dyu{s} = ErrorEstimates.ZMAGCL95u(:,comp) - y{s};
            end
            if popts.type == 2
                f = S{s}.fe*S{s}.Metadata.freqsf;
                y{s} = z2rho(f, S{s}.Z(:,comp));
                dyl{s} = y{s} - ErrorEstimates.RHOCL95l(:,comp);
                dyu{s} = ErrorEstimates.RHOCL95u(:,comp) - y{s};
            end
            if popts.type == 3
                y{s} = real(S{s}.Z(:,comp));
                dyl{s} = y{s} - real(ErrorEstimates.ZCL95l(:,comp));
                dyu{s} = real(ErrorEstimates.ZCL95u(:,comp)) - y{s};
            end
        end

        if strcmp(panel,'bottom')
            if popts.type == 3
                y{s} = imag(S{s}.Z(:,comp));
                dyl{s} = y{s} - imag(ErrorEstimates.ZCL95l(:,comp));
                dyu{s} = imag(ErrorEstimates.ZCL95u(:,comp)) - y{s};
            else
                ang = atan2(imag(S{s}.Z(:,comp)),real(S{s}.Z(:,comp)));
                if popts.unwrap
                    y{s} = (180/pi)*unwrap(ang);
                else
                    y{s} = (180/pi)*ang;
                end
                dyl{s} = y{s} - ErrorEstimates.PHICL95l(:,comp);
                dyu{s} = ErrorEstimates.PHICL95u(:,comp)- y{s};
            end
        end

        if popts.vs_period
            x{s} = 1./(S{s}.fe*S{s}.Metadata.freqsf);
        else
            x{s} = S{s}.fe*S{s}.Metadata.freqsf;
        end

        if popts.vs_period && ~isempty(popts.period_range)
            idx = x{s} <= popts.period_range(1) | x{s} >= popts.period_range(2);
            y{s}(idx) = NaN;
        end
    end

    if length(S) == 1
        x = x{s};
        y = y{s};
        dyl = dyl{s};
        dyu = dyu{s};
    end
end % function xyvals_()

function str = unitstr_(meta)

    str = '';
    if isempty(meta.inunit) || isempty(meta.outunit)
        return
    end
    inunit = sprintf('%s',meta.inunit);
    if contains(inunit,'/')
        inunit = sprintf('(%s)',inunit);
    end
    outunit = sprintf('%s',meta.outunit);
    if contains(outunit,'/')
        outunit = sprintf('(%s)',outunit);
    end
    str = sprintf('[%s/%s]',outunit,inunit);
end % function unitstr_()