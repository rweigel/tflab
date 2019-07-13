function adjust_ylim()

yt = get(gca,'YTick');
yl = get(gca, 'YLim');        
% Expand axis limits so legend does not overlap data.
yl(end) = yl(end) + (yt(2)-yt(1))/2;
if (yl(1) == -yl(end))
    % If limits were symmetric, adjust bottom to keep symmetry.
    yl(1) = yl(1) - (yt(2)-yt(1))/2;
end
set(gca, 'YLim', yl);

