function adjust_ylim(pos)
% ADJUST_YLIM
%
%   Expand axis limits so horizontal legend or vertical line labels are
%   less likely to overlap data.

if nargin == 0
    pos = 'upper';
end

drawnow;
yt = get(gca,'YTick');
yl = get(gca,'YLim');

if strcmp(pos, 'upper') || strcmp(pos, 'both')
    if strcmp(get(gca(),'YScale'),'log')
        if length(yt) > 1
            yl(end) = 10^(log10(yl(end)) + 2*(log10(yt(end))-log10(yt(end-1))));
        end
    else
        if (yl(1) == -yl(end)) && ~strcmp(pos, 'both')
            % If limits were symmetric, adjust lower to keep symmetry.
            yl(1) = yl(1) - (yt(end)-yt(end-1));
        end
        yl(end) = yl(end) + (yt(end)-yt(end-1));
    end
end
if strcmp(pos,'lower') || strcmp(pos, 'both')
    if strcmp(get(gca(),'YScale'),'log')
        if length(yt) > 1
            yl(1) = 10^(log10(yl(1)) - 0.5*(log10(yt(2))-log10(yt(1))));
        end
    else
        if (yl(1) == -yl(end)) && ~strcmp(pos, 'both')
            % If limits were symmetric, adjust upper to keep symmetry.
            yl(end) = yl(end) - (yt(2)-yt(1));
        end
        yl(1) = yl(1) - (yt(2)-yt(1));
    end
end
set(gca,'YLim',yl);
set(gca,'YTick',yt);
drawnow;
