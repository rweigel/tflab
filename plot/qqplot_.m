function qqplot_(S,fidx,cidx,sidx)
%QQPLOT_ Quantitle-Quantile plot
%
%   QQPLOT_(S)
%   QQPLOT_(S,freq_idx)
%   QQPLOT_(S,freq_idx,comp_idx)
%   QQPLOT_(S,freq_idx,comp_idx,seg_idx)
%
%   where freq_idx is index of S.fe, comp_idx is column of Z, and
%   seg_idx, is the segment number (if applicable).

opts = struct();
    opts.title = '';
    opts.print = 0;
    opts.printname = 'qqplot';
    opts.printdir = '';
    opts.printfmt = {'pdf'};

if nargin == 1
    fidx = 3;
end
if nargin < 3
    cidx = 1;
end
if nargin < 4
    sidx = -1;
end

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S)
    fe = S{1}.fe;
    for i = 1:length(S)
        if length(fe) ~= S{i}.fe | all(fe == S{i}.fe)
            
        end
    end
    fe = fe(fidx);
else
    fe = S.fe(fidx);
end

legendstr = {};
if iscell(S)
    for i = 1:length(S)
        if isfield(S{i},'Regression')
            Residuals{i} = S{i}.Regression.Residuals{fidx,cidx};
        else
            if sidx == -1
                Residuals{i} = cat(2,S{i}.Segment.Regression.Residuals{fidx,cidx,:});
            else
                Residuals{i} = S{i}.Segment.Regression.Residuals{fidx,cidx,sidx};
            end
        end
        legendstr{i} = S{i}.Options.description;            
    end
    titlestr = sprintf('$f=$ %g; $T=$ %.3f',fe,1/fe);
else
    if isfield(S,'Regression')
        Residuals = S.Regression.Residuals{fidx,cidx};
    else
        if sidx == -1
            Residuals = cat(2,S.Segment.Regression.Residuals{fidx,cidx,:});
        else
            Residuals = S.Segment.Regression.Residuals{fidx,cidx,sidx};
        end
    end
    titlestr = sprintf('%s\n$f=$ %g; $T=$ %.3f',S.Options.description,fe,1/fe);
end

% TODO: Code copied from zplot.m
if (iscell(S) && size(S{1}.Z,2) > 1) || (isstruct(S) && size(S.Z,2) > 1)
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{yx}','\rho^a_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Rhostrs = {'\rho^a'};
    Phistrs = {'\phi'};
end

if ~isfield(S,'Regression') && sidx ~= -1
    norm = std(S.Segment.Regression.Residuals{fidx,cidx,sidx});
    titlestr = sprintf('%s\nSeg. = %d/%d; SN = %.1f; Coh = %.2f; $\\sigma_{\\Delta %s}$ = %.1e',...
        titlestr,sidx,size(S.Segment.Regression.Residuals,3),...
        S.Segment.Metrics.SN.Smoothed(fidx,cidx,sidx),...
        S.Segment.Metrics.Coherence.Smoothed(fidx,cidx,sidx),...
        Zstrs{cidx},norm);
end

PositionTop = [0.1300 0.5400 0.7750 0.38];
PositionBottom = [0.1300 0.1100 0.7750 0.38];

figprep();
ax1 = subplot('Position',PositionTop);
    grid on;box on;hold on;axis square;
    plotqq(Residuals,'real')
    hold on;grid on;box on;
    title(titlestr);
    xlabel('');
    set(gca,'XTickLabels','');
    ylabel(sprintf('Re[$\\Delta %s$] Quantiles',Zstrs{cidx}));
    if ~isempty(legendstr)
        legend(legendstr,'Location','SouthEast');
    end
    adjust_exponent();    
ax2 = subplot('Position',PositionBottom);
    grid on;box on;hold on;axis square;
    plotqq(Residuals,'imag');
    hold on;grid on;box on;    
    xlabel('Standard Normal Quantiles');
    ylabel(sprintf('Im[$\\Delta %s$] Quantiles',Zstrs{cidx}));
    adjust_exponent()
    legend off;

%return
% Plot diagonal line    
axes(ax1);hold on;
xl1 = get(ax1,'XLim');    
yl1 = get(ax1,'YLim');
axes(ax2);hold on;
xl2 = get(ax2,'XLim');
yl2 = get(ax2,'YLim');
m = max(abs([xl1,yl1,xl2,yl2]));
plot([-m,m],[-m,m],'k');    
axes(ax1)
plot([-m,m],[-m,m],'k');    

function plotqq(Residuals,comp)
    if strcmp(comp,'real')
        compdata = @(x) real(x);
    end
    if strcmp(comp,'imag')
        compdata = @(x) imag(x);
    end
    if iscell(Residuals)
        colors = [0,0,0; lines(length(Residuals))];
        for i = 1:length(Residuals)
            % Hack to get correct legend symbol colors.
            plot(NaN, NaN,'Color',colors(i,:),...
                'LineStyle','none','Marker','.','MarkerSize', 10);
        end
        for i = 1:length(Residuals)
            %qqplot(compdata(Residuals{i}));
            [x,y] = qqdata(compdata(Residuals{i}));
            plot(x,y,'Color',colors(i,:),...
                'LineStyle','none','Marker','.','MarkerSize', 10);
        end
    else
        %qqplot(compdata(Residuals));
        [x,y] = qqdata(compdata(Residuals));
        plot(x,y,'k.','MarkerSize', 10);
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
    x = 1:size(yo,1);
    x = (x-.5)./nvalid;
    x(isnan(yos)) = NaN;  
    x = norminv(x);
end
end