function adjust_ylim(pos)
% ADJUST_YLIM
%
%   Expand axis limits so horizontal legend or vertical line labels are
%   less likely to overlap data.
%
%   ADJUST_YLIM() Adjust upper and lower limits
%   ADJUST_YLIM('both') Adjust upper and lower limits
%   ADJUST_YLIM('lower') Only adjust lower limit
%   ADJUST_YLIM('upper) Only adjust upper limit
%

% TODO: One should be able to compute exact adjustment needed.

direction = 'y';

if nargin == 0
    pos = 'upper';
end

debug = 0;

ax = gca();
drawnow;

ytick = get(ax,'YTick');
ytick_old = ytick;

ylim = get(ax,'YLim');

if ~isempty(ax.UserData) && isfield(ax.UserData, 'YLimLast')
    tf = all(ax.UserData.YLimLast == ylim);
    if all(tf)
        if debug
            fprintf('adjust_ylim(): Previous limits are same as current.\n');
        end
        return;
    end
end

if debug
    fprintf('adjust_ylim():\n');
    ytick
    ylim
end

if strcmp(pos, 'upper') || strcmp(pos, 'both')
    if strcmp(get(ax,'YScale'),'log')
        if length(ytick) > 1
            ylim(end) = ylim(end)*10^(0.5*(log10(ytick(end))-log10(ytick(end-1))));
        else
            ylim(end) = 10*ytick;
        end
    else
        if (ylim(1) == -ylim(end)) && ~strcmp(pos, 'both')
            % If limits were symmetric, adjust lower to keep symmetry.
            ylim(1) = ylim(1) - (ytick(end)-ytick(end-1));
        end
        ylim(end) = ylim(end) + (ytick(end)-ytick(end-1));
    end
end

if strcmp(pos,'lower') || strcmp(pos, 'both')
    if strcmp(get(ax,'YScale'),'log')
        if length(ytick) > 1
            ylim(1) = ylim(1)/10^(0.5*(log10(ytick(2))-log10(ytick(1))));
            ylim(1) = 0.5*10^(log10(ylim(1)) - 0.5*(log10(ytick(2))-log10(ytick(1))));
        else
           ylim(end) = ytick/10'
        end
    else
        if (ylim(1) == -ylim(end)) && ~strcmp(pos, 'both')
            % If limits were symmetric, adjust upper to keep symmetry.
            ylim(end) = ylim(end) - (ytick(2)-ytick(1));
        end
        ylim(1) = ylim(1) - (ytick(2)-ytick(1));
    end
end

if ylim(1) > ylim(2)
    return;
end

if debug
    fprintf('adjust_ylim(): Setting ax.UserData.YLimIgnoreChange = 1\n');
end
ax.UserData.YLimIgnoreChange = 1;

if debug
    ylim
end

set(ax,'YLim',ylim);
if length(ytick) > 1
    set(ax,'YTick',ytick_old);
end
drawnow;

if debug
    fprintf('adjust_ylim(): Setting ax.UserData.YLimIgnoreChange = 0\n');
end

ax.UserData.YLimIgnoreChange = 0;
ax.UserData.YLimLast = ylim;

if isprop(ax.YAxis,'LimitsChangedFcn')
    ax.YAxis.LimitsChangedFcn = @(src,evt) reset(src,evt,debug);
else
    addlistener(gca(), [upper(direction),'Lim'], 'PostSet', @(obj,evt) reset(obj,evt,debug));
end

function reset(obj,evt,debug)
    if debug
        fprintf(['adjust_ylim(): Reset called. Setting TickLabelMode to auto '...
                 'and deleting listener for %sLim change.\n'],direction);
    end
    set(gca,[upper(direction), 'TickMode'],'auto');
    delete(obj);
end
end
