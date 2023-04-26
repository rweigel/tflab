function adjust_yticks(t,ax,no_lead_number)
% ADJUST_YTICKS - Change labels of yticks when increment is small

direction = 'y';

debug = 0;

if nargin == 0
    t = 1e-4; % Threshold to apply adjustment.
end

if nargin < 2
    ax = gca;
else
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
end

if nargin < 3
    no_lead_number = 1;
end

yt = get(ax,'YTick')';
dy = diff(yt);

if min(dy) > t
    return
end
if strcmp(get(gca,'YScale'),'log')
    %set(gca,'YScale','linear');
    return;
    if ~all(dy < t)
        %return
    end
end

zero_in_interval = 0;
if max(yt) > 0 && min(yt) < 0
    zero_in_interval = 1;
end

drawnow;
yl = get(ax,'YTickLabels');

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
        fprintf('Original: %s\n',yl{i});
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
        fprintf('New:      %s\n',replace(replace(yl{i},hspace,''),'\cdot',''));
    end
end

set(gca,'TickLabelInterpreter','latex');
set(gca,'YTickLabels',yl);

ax = gca();

if isprop(ax.YAxis,'LimitsChangedFcn')
    ax.YAxis.LimitsChangedFcn = @(src,evt) reset(src,evt,debug);
else
    addlistener(gca(), [upper(direction),'Lim'], 'PostSet', @(obj,evt) reset(obj,evt,debug));
end

function reset(obj,evt,debug)
    if debug
        fprintf(['Reset called. Setting TickLabelMode to auto '...
                 'and deleting listener for %sLim change.\n'],direction);
    end
    set(gca,[upper(direction), 'TickLabelMode'],'auto');
    delete(obj);
    adjust_yticks(t,ax,no_lead_number);
end
end

