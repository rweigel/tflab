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

S = tflab_metadata(S);

% TODO: Check all same units and sfs or allow different.
frequnit = S{1}.Metadata.frequnit;
popts   = tflabplot_options(S, popts, 'zplot');

if nargin < 3
    for s = 1:length(S)
        nr(s) = size(popts.zstrs{s},1);
        nc(s) = size(popts.zstrs{s},2);
        logmsg('TF #%d is %dx%d\n',s,nr(s),nc(s));
    end
    % Compute list of all possible components to plot
    cn = 1;
    for i = 1:max(nr)
        for j = 1:max(nc)
            comps{cn} = [i,j];
            cn = cn + 1;
        end
    end
end

% Single transfer function
if length(S) == 1
    logmsg('Plotting single transfer function.\n')

    yl1 = unitstr_(S{1}.Metadata);
    yl2 = '[$^\circ$]';
    if popts.type == 1 || popts.type == 2
        template2 = '$%s$';
        strings2 = popts.phistrs{1};
        if popts.type == 1
            template1 = '$|%s|$';
            strings1 = popts.zstrs{1};
        else
            template1 = '$%s$';
            strings1 = popts.rhostrs{1};
        end
    else
        template1 = 'Re$(%s)$';
        strings1 = popts.zstrs{1};
        template2 = 'Im$(%s)$';
        strings2 = popts.zstrs{1};
        yl2 = unitstr_(S{1}.Metadata);
    end

    for c = 1:length(comps)
        idx = sub2ind(size(popts.zstrs{1}), comps{c}(1), comps{c}(2));
        ls1{c} = sprintf(template1,strings1{idx})
        [x1,y1(:,c),dyu1(:,c),dyl1(:,c)] = xyvals_(S,popts,[comps{c}(1), comps{c}(2)],'top','parametric');
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            [~,~,dyu1b(:,c),dyl1b(:,c)] = xyvals_(S,popts,[comps{c}(1), comps{c}(2)],'top','bootstrap');
        end

        ls2{c} = sprintf(template2,strings2{idx});
        [x2,y2(:,c),dyu2(:,c),dyl2(:,c)] = xyvals_(S,popts,[comps{c}(1), comps{c}(2)],'bottom','parametric');
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            [~,~,dyu2b(:,c),dyl2b(:,c)] = xyvals_(S,popts,[comps{c}(1), comps{c}(2)],'bottom','bootstrap');
        end
    end

    figprep();
    ax1 = subplot('Position', popts.PositionTop);

        popts.line = markeropts(size(S{1}.Z,1), 1);

        h1 = plot(x1, y1, popts.line{:});
        hold on;grid on;box on;
        title_(S{1},popts,'z');
        colororder_(ax1,y1);

        if length(ls1) > 1
            ylabel(yl1);
            %legend(ls1,popts.legend{:});
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
            errorbars_(h1,0.97*x1,y1,dyl1,dyu1,1);
            errorbars_(h1,x1*1.03,y1,dyl1b,dyu1b,1);
        else
            errorbars_(h1,x1,y1,dyl1,dyu1,2);
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

        if popts.type ~= 3 && ~popts.unwrap
            set(gca,'YScale','linear');
            set(gca,'YLim',[-180,180]);
            set(gca,'YTick',-180:45:180);
        end
        if popts.type == 1
            adjust_ylim('both');
        end

        % Must be called after scale type is set.
        if isfield_(S{1},'Metrics.ErrorEstimates.Bootstrap')
            errorbars_(h2,0.97*x2,y2,dyl2,dyu2,1);
            errorbars_(h2,1.03*x2,y2,dyl2b,dyu2b,1);
        else
            errorbars_(h2,x2,y2,dyl2,dyu2,2);
        end

        %adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,1,frequnit);

        xlims = get(gca(),'XLim');
        line(xlims,[0,0])

    if popts.print == 1
        ext = '';
        for c = 1:length(comps)
           idx = sub2ind(size(popts.zstrs{1}), comps{c}(1), comps{c}(2));
           ext = [ext,popts.zstrs{1}{idx}];
        end
        figsave_(popts,ext);
    end

end % if isstruct(S)

% Multiple transfer functions
if length(S) > 1

    if length(comps) == 1
        comp = comps{1};
        logmsg('Plotting component Z_%d%d for TFs with this component.\n',comp(1),comp(2));
    else
        msg = 'Plotting %d components for %d transfer functions.\n';
        logmsg(msg,length(comps),length(S));
        for c = 1:length(comps)
            if c > 1
                figure;
            end
            [ax1(c),ax2(c)] = zplot(S,popts,{[comps{c}(1), comps{c}(2)]});
        end
        return
    end

    figprep();
    ax1 = subplot('Position', popts.PositionTop);
        [x,y,dyu,dyl,kept] = xyvals_(S,popts,comp,'top','parametric');

        for k = 1:length(kept)
            s = kept(k);
            ls{k} = sprintf('%s',S{s}.Options.description);
            popts.line = markeropts(size(S{s}.Z,1), s);
            h(k) = plot(x{s}, y{s}, popts.line{:});
            if k == 1
                grid on;hold on;box on;
            end
        end
        idx = sub2ind(size(popts.zstrs{s}), comp(1), comp(2));
        legend(h,ls,popts.legend{:});
        title_(S,popts,'z');

        if popts.type < 3
            set(gca,'YScale','log');
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        yunitstr_ = unitstr_(S{1}.Metadata);
        if popts.type == 1
            yl = sprintf('$|%s|$%s',popts.zstrs{s}{idx},yunitstr_);
            ylraw = popts.zstrs{s}{idx};
        end
        if popts.type == 2
            yl = sprintf('$%s$ [$\\Omega\\cdot$m]',popts.rho{s}{idx});
            ylraw = popts.rho{s}{idx};
        end
        if popts.type == 3
            yl = sprintf('Re$(%s)$ %s',popts.zstrs{s}{idx},yunitstr_);
            ylraw = popts.zstrs{s}{idx};
        end
        if ~isempty(popts.period_range)
            set(gca,'XLim',popts.period_range);
        end
        ylabel(yl);

        % ebars_() must be called after scale type is set.
        ebars_(h, x(kept), y(kept), dyl(kept), dyu(kept), 1, ylraw)

        adjust_yticks(1e-4);
        adjust_exponent('y');

        if popts.type == 3
            adjust_ylim('both');
        else
            adjust_ylim('upper');
        end
        setx(popts,0,frequnit);

    ax2 = subplot('Position', popts.PositionBottom);

        [x,y,dyu,dyl,kept] = xyvals_(S,popts,comp,'bottom','parametric');

        for k = 1:length(kept)
            s = kept(k);
            popts.line = markeropts(size(S{s}.Z,1), s);
            if popts.vs_period
                h(k) = semilogx(x{s}, y{s}, popts.line{:});
            else
                h(k) = plot(x{s}, y{s}, popts.line{:});
            end

            if k == 1,grid on;box on;hold on;end
        end
        idx = sub2ind(size(popts.zstrs{s}), comp(1), comp(2));
        yl = sprintf('$%s$ [$^\\circ]$',popts.phistrs{s}{idx});
        ylraw = popts.phistrs{s}{idx};
        if popts.type == 3
            yl = sprintf('Im$(%s)$ %s',popts.zstrs{s}{idx},unitstr_(S{1}.Metadata));
            ylraw = popts.zstrs{s}{idx};
        end
        ylabel(yl);
        
        if popts.type ~= 3 && ~popts.unwrap
            set(gca,'YScale','linear');
            set(gca,'YLim',[-180,180]);
            set(gca,'YTick',-180:45:180);
        end
        if popts.type == 1
            adjust_ylim('both');
        end
        
        setx(popts,1,frequnit);
        ebars_(h, x(kept), y(kept), dyl(kept), dyu(kept), 1, ylraw)

        adjust_yticks(1e-4);
        adjust_exponent();

    var_name = replace(popts.zstrs{s}{idx},'\delta ','delta');
    figsave_(popts,var_name);

end % if iscell(S)

end % function zplot()

function ebars_(h, x, y, dyl, dyu, lw, yl)
    if isscalar(h)
        errorbars_(h(1),x{1},y{1},dyl{1},dyu{1},lw);
    elseif length(h) < 5
        %errorbars_(h,{0.97*x{1},1.03*x{2}},{y{1},y{2}},{dyl{1},dyl{2}},{dyu{1},dyu{2}},lw);
        %errorbars_(h(2),1.03*x{2},y{2},dyl{2},dyu{2},lw);
        errorbars_(h, x, y, dyl, dyu, lw);
    else
        logmsg('Not plotting error bars for %s because more than four TFs being plotted.\n', yl);
    end
end

function [x,y,dyu,dyl,kept] = xyvals_(S,popts,comp,panel,errorbar_method)

    % TODO: Assumes units are the same for all Zs.
    %       Check this before plotting.

    kept = [];

    for s = 1:length(S)
        if comp(1) > size(popts.zstrs{s}, 1) || comp(2) > size(popts.zstrs{s}, 2)
            continue
        end
        kept(end+1) = s;
        % Matrixes have columns of Z11, Z12, Z13, ..., Z21, Z22, Z23, ...
        nr = size(popts.zstrs{s}, 1); % number of rows
        nc = size(popts.zstrs{s}, 2); % number of columns
        idx = (comp(1) - 1)*nc + comp(2);
        if strcmp(errorbar_method,'parametric')
            ErrorEstimates = S{s}.Metrics.ErrorEstimates.Parametric;
        end
        if strcmp(errorbar_method,'bootstrap')
            ErrorEstimates = S{s}.Metrics.ErrorEstimates.Bootstrap;
        end
        if strcmp(panel,'top')
            if popts.type == 1
                y{s} = abs(S{s}.Z(:,idx));
                dyl{s} = y{s} - ErrorEstimates.ZMAGCL95l(:,idx);
                dyu{s} = ErrorEstimates.ZMAGCL95u(:,idx) - y{s};
            end
            if popts.type == 2
                f = S{s}.fe*S{s}.Metadata.freqsf;
                y{s} = z2rho(f, S{s}.Z(:,idx));
                dyl{s} = y{s} - ErrorEstimates.RHOCL95l(:,idx);
                dyu{s} = ErrorEstimates.RHOCL95u(:,idx) - y{s};
            end
            if popts.type == 3
                y{s} = real(S{s}.Z(:,idx));
                dyl{s} = y{s} - real(ErrorEstimates.ZCL95l(:,idx));
                dyu{s} = real(ErrorEstimates.ZCL95u(:,idx)) - y{s};
            end
        end

        if strcmp(panel,'bottom')
            if popts.type == 3
                y{s} = imag(S{s}.Z(:,idx));
                dyl{s} = y{s} - imag(ErrorEstimates.ZCL95l(:,idx));
                dyu{s} = imag(ErrorEstimates.ZCL95u(:,idx)) - y{s};
            else
                ang = atan2(imag(S{s}.Z(:,idx)),real(S{s}.Z(:,idx)));
                if popts.unwrap
                    y{s} = (180/pi)*unwrap(ang);
                else
                    y{s} = (180/pi)*ang;
                end
                dyl{s} = abs(y{s} - ErrorEstimates.PHICL95l(:,idx));
                dyu{s} = dyl{s};
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