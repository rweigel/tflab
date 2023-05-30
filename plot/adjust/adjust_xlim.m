function adjust_xlim()
%ADJUST_XLIM
%
% Expand axis limits so small gap between bounding box and first tick.
% (Use when data have spike at start or end that would be difficult to see
% otherwise.)

xt = get(gca,'XTick');
xl = get(gca,'XLim');        

xl(end) = xl(end) + (xt(2)-xt(1))/10;
xl(1) = xl(1) - (xt(2)-xt(1))/10;

set(gca, 'XLim', xl);