function [ax1, ax2] = qqplot_(S,popts,cidxs,fidxs)
%QQPLOT_  Quantitle-Quantile and residuals plots
%
%   QQPLOT_(S)
%   QQPLOT_(S, opts)
%
%   If popts.type is 'standard' a single plot is created for each component and
%   frequency. If popts.type is a single "rotated QQ" plot is created for each
%   component and all frequencies.
%
%   QQPLOT_(S,opts,comp_idx)
%   QQPLOT_(S,opts,comp_idx,freq_idx)
%
%   where comp_idx is column component of output and freq_idx is the index of fe.
%

msg = 'S must be a tflab struct or cell array of tflab structs';
assert(isstruct(S) || iscell(S), msg);

if nargin < 2
    popts = struct();
end
popts = tflabplot_options(S, popts, 'qqplot');

if ~iscell(S)
    S = {S};
end
if nargin < 5
    sidx = -1;
end

if nargin < 3
    % Assumes residuals have the same number of components for all freqs.
    cidxs = 1:size(S{1}.Metrics.Residuals{1},2);
end
if nargin < 4
    fidxs = 1:size(S{1}.Metrics.Residuals,1);
end

logmsg('Plotting\n');

if strcmp(popts.type,'standard') && length(fidxs) > 1 && length(cidxs) > 1
    first = 1;
    for cidx = cidxs
        for fidx = fidxs
            if first == 0
                figure();figprep();
            end
            first = 0;
            qqplot_(S{1},popts,cidx,fidx);
        end
    end
    return
end

frequnit = S{1}.Metadata.frequnit;
timeunit = S{1}.Metadata.timeunit;

if strcmp(popts.type,'standard')
    fidx = fidxs(1);
    cidx = cidxs(1);
    outstr = replace(popts.outstr{cidx},'$','');

    legendstr = {};
    for i = 1:length(S)
        Residuals{i} = S{i}.Metrics.Residuals{fidx}(:,cidx);
        fe{i} = S{i}.Metrics.fe(fidx,1);
        legendstr{i} = S{i}.Options.description;
    end
    % TODO: Check that all fes are the same.
    titlestr = sprintf('$f=$ %g [%s]; $T=$ %.3f [%s]',...
                        fe{1},frequnit,1/fe{1},timeunit);

    figprep();
    ax1 = subplot('Position',popts.PositionTop);
        grid on;box on;hold on;axis square;
        plotqq(Residuals,'real');
        colororder_(ax1,Residuals);
        hold on;grid on;box on;
        title(titlestr);
        xlabel('');
        set(gca,'XTickLabels','');
        ylabel(sprintf('Re[$\\widetilde{\\Delta %s}$] Quantiles',outstr));
        if ~isempty(legendstr)
            legend(legendstr,'Location','NorthWest','box','off','color','none');
        end
        adjust_exponent();
    ax2 = subplot('Position',popts.PositionBottom);
        grid on;box on;hold on;axis square;
        plotqq(Residuals,'imag');
        colororder_(ax2,Residuals);
        hold on;grid on;box on;
        xlabel('Standard Normal Quantiles');
        ylabel(sprintf('Im[$\\widetilde{\\Delta %s}$] Quantiles',outstr));
        adjust_exponent()
        legend off;

    % Plot diagonal line
    xl1 = get(ax1,'XLim');
    yl1 = get(ax1,'YLim');
    xl2 = get(ax2,'XLim');
    yl2 = get(ax2,'YLim');
    m = ceil(max(abs([xl1,yl1,xl2,yl2])));
    axes(ax1)
        plot([-m,m],[-m,m],'k');
    axes(ax2)
        plot([-m,m],[-m,m],'k');

    set([ax1,ax1],'XLim',[-m,m],'YLim',[-m,m]);
    set([ax1,ax1],'YTick',[-m:m]);
    set([ax1,ax1],'XTick',[-m:m]);
    set([ax2,ax2],'XLim',[-m,m],'YLim',[-m,m]);
    set([ax2,ax2],'YTick',[-m:m]);
    set([ax2,ax2],'XTick',[-m:m]);

    ext = sprintf('%s-fidx_%d',outstr,fidx);
    figsave_(popts,ext)
end

if strcmp(popts.type,'combined')
    if length(cidxs) > 1
        for cidx = cidxs
            if cidx > 1
                figure();figprep();
            end
            qqplot_(S,popts,cidx);
        end
        return
    else
        cidx = cidxs(1);
    end
    outstr = replace(popts.outstr{cidx},'$','');
    for fidx = 1:size(S{1}.Metrics.Residuals,1) % frequencies
        Residuals = S{1}.Metrics.Residuals{fidx}(:,cidx);
        fe{fidx} = S{1}.Metrics.fe(fidx,1);
        [xr{fidx},yr{fidx}] = qqdata_rotated(Residuals,'real');
        [xi{fidx},yi{fidx}] = qqdata_rotated(Residuals,'imag');
    end
    figprep();
    popts.vs_period = 1;
    ax1 = subplot('Position',popts.PositionTop);
        for fidx = 1:length(fe)
            semilogx(10.^xr{fidx}+1./fe{fidx},yr{fidx},'k.','MarkerSize', 10);
            hold on;
        end
        ylabel(sprintf('Re[$\\widetilde{\\Delta %s}$] QQ parallel',outstr));
        setx(popts,0,frequnit);
    ax2 = subplot('Position',popts.PositionBottom);
        for fidx = 1:length(fe)
            semilogx(10.^xi{fidx}+1./fe{fidx},yi{fidx},'k.','MarkerSize', 10);
            hold on;
        end
        ylabel(sprintf('Im[$\\widetilde{\\Delta %s}$] QQ parallel',outstr));
        setx(popts,1,frequnit);
        xl = get(gca,'XLabel');
        txt = sprintf('%s and %s',xl.String,'QQ perpendicular');
        ax2.XLabel.String = txt;
    ext = sprintf('%s-combined',outstr);
    figsave_(popts,ext)
end


if 0 && ~isfield(S,'Regression') && sidx ~= -1
    norm = std(S.Segment.Regression.Residuals{fidx,cidx,sidx});
    titlestr = sprintf('%s\nSeg. = %d/%d; SN = %.1f; Coh = %.2f; $\\sigma_{\\Delta %s}$ = %.1e',...
        titlestr,sidx,size(S.Segment.Regression.Residuals,3),...
        S.Segment.Metrics.SN(fidx,cidx,sidx),...
        S.Segment.Metrics.Coherence(fidx,cidx,sidx),...
        Zstrs{cidx},norm);
end

end

function plotqq(Residuals,comp)
    if strcmp(comp,'real')
        compdata = @(x) real(x);
    end
    if strcmp(comp,'imag')
        compdata = @(x) imag(x);
    end
    if iscell(Residuals)
        for i = 1:length(Residuals)
            % Hack to get correct legend symbol colors.
            %plot(NaN, NaN,'LineStyle','none','Marker','.','MarkerSize', 10);
        end
        for i = 1:length(Residuals)
            %qqplot(compdata(Residuals{i}));
            [x,y] = qqdata(compdata(Residuals{i}));
            plot(x,y,'LineStyle','none','Marker','.','MarkerSize', 10);
        end
    else
        %qqplot(compdata(Residuals));
        [x,y] = qqdata(compdata(Residuals));
    end
end

function [x,y] = qqdata(yo)

    % TODO: Could vectorize
    if size(yo,2) > 1
        for j = 1:size(yo,2)
            [x(:,j),y(:,j)] = qqdata(yo(:,j));
        end
        return
    end

    % Based on qqplot.m
    yos = sort(yo);
    y = yos/std(yos(~isnan(yo)));

    nvalid = sum(~isnan(yos));
    x = [1:size(yo,1)]';
    x = (x-.5)./nvalid;
    x(isnan(yos)) = NaN;
    x = norminv(x);
end

function [x,y] = qqdata_rotated(residuals,comp)
    if strcmp(comp,'real')
        [x,y] = qqdata(real(residuals));
    else
        [x,y] = qqdata(imag(residuals));
    end
    % Rotate data 45 degrees ccw
    r = hypot(x,y);
    theta = atan2(y,x);
    x = -r.*sin(theta-pi/4);
    y = r.*cos(theta-pi/4);
end
