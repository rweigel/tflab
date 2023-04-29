function adjust_yticks(t,ax,no_lead_number, listen)
% ADJUST_YTICKS - Change labels of yticks when increment is small

return;

direction = 'y';

debug = 1;

if nargin < 4
    listen = 1;
end

if nargin == 0
    t = 1e-4; % Threshold to apply adjustment.
end

if nargin < 3
    no_lead_number = 1;
end

ax = gca();

set(ax, 'YTickMode', 'auto', 'YTickLabelMode', 'auto');

drawnow;
yl = get(ax,'YTickLabels');
yt = get(ax,'YTick')';
dy = diff(yt);

if debug
    yt
    yl
end

if ~isempty(ax.UserData) && isfield(ax.UserData,'YTickLabelsLast')
    if length(ax.UserData.YTickLabelsLast) == length(yl)
        tf = strcmp(ax.UserData.YTickLabelsLast, yl);
        tf
        if all(tf)
            if debug
                fprintf('adjust_yticks(): Previous labels are same as current.\n');
            end
            return;
        end
    end
end

if min(dy) > t
    return;
end

if strcmp(get(gca(),'YScale'),'log')
    set(gca(),'YScale','linear');
    drawnow;
    yl = get(ax,'YTickLabels');
    yt = get(ax,'YTick')';    
end

zero_in_interval = 0;
if max(yt) > 0 && min(yt) < 0
    zero_in_interval = 1;
end

% Prevent number of ticks from changing if vertical window size changes.
% There is no callback for this, and if it happens, labels may become
% incorrect.
yticks('manual');

if strcmp(yl{1}(1),'$')
    % Already modified
    %return
end

for i = 1:length(yl)
    lens(i) = length(yl{i});
end

[~,idx] = min(lens);
dym = min(dy);
ed = floor(log10(dym)); % Exponent digit

hspace = '\hspace{-0.25em}';

for i = 1:length(yl)
    if i == idx
        continue
    end
    if yt(i)-yt(idx) > 0
        sign = sprintf('%s + %s',hspace,hspace);
        dy = floor((yt(i)-yt(idx))/dym);
    else
        sign = sprintf('%s - %s',hspace,hspace);
        dy = floor((yt(idx)-yt(i))/dym);  
    end
    dy = dy*dym/10^ed;
    if debug
        fprintf('Orig: %s\n',yl{i});
    end
    if dy == 1
        yl{i} = sprintf('$%s%s10^{%d}$',yl{idx},sign,ed);
    else
        if no_lead_number || zero_in_interval
            % \phantom{} is needed so first hspace is evaluated.
            yl{i} = sprintf('$\\phantom{}%s%.0f%s\\cdot%s 10^{%d}$',...
                            sign,dy,hspace,hspace,ed);
        else
            yl{i} = sprintf('$%s%s%.0f%s\\cdot%s 10^{%d}$',...
                            yl{idx},sign,dy,hspace,hspace,ed);
        end                        
    end
    if debug
        fprintf('New:  %s\n',replace(replace(yl{i},hspace,''),'\cdot',''));
    end
end

set(ax,'TickLabelInterpreter','latex');
set(ax,'YTickLabels',yl);
drawnow;
ax.UserData.YTickLabelsLast = yl;

if listen == 0
    return;
end

if 0 && isprop(ax.YAxis,'LimitsChangedFcn')
    ax.YAxis.LimitsChangedFcn = @(src,evt) reset(src,evt,debug);
else
    addlistener(gca(), [upper(direction),'Lim'], 'PostSet', @(obj,evt) reset(obj,evt,debug));
end

function reset(obj,evt,debug)
    if debug
        fprintf(['adjust_yticks(): Reset called. Setting TickLabelMode to auto '...
                 'and deleting listener for %sLim change.\n'],direction);
    end
    set(ax,[upper(direction), 'TickLabelMode'],'auto');
    delete(obj);
    adjust_yticks(t,ax,no_lead_number,0);
end
end

