function adjust_ylim(pos)
% ADJUST_YLIM
%
%   Expand axis limits so horizontal legend or vertical line labels are
%   less likely to overlap data.

if nargin == 0
    pos = 'upper';
end

if strcmp(pos, 'both')
    adjust_ylim('upper');
    adjust_ylim('lower');
    return
end

yt = get(gca,'YTick');
yl = get(gca,'YLim');

if strcmp(pos, 'upper')
    if strcmp(get(gca(),'YScale'),'log')
        yl(end) = 10*yl(end);
    else
        yl(end) = yl(end) + (yt(end)-yt(end-1))/2;
        if (yl(1) == -yl(end))
            % If limits were symmetric, adjust lower to keep symmetry.
            yl(1) = yl(1) - (yt(end)-yt(end-1))/2;
        end
    end
else
    if strcmp(get(gca(),'YScale'),'log')
        yl(1) = yl(1)/10;
    else
        yl(1) = yl(1) - (yt(2)-yt(1))/2;
        if (yl(1) == -yl(end))
            % If limits were symmetric, adjust upper to keep symmetry.
            yl(end) = yl(end) - (yt(2)-yt(1))/2;
        end
    end
end
set(gca,'YLim',yl);
