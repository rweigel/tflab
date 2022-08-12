function adjust_yticks(t)


if nargin == 0
    t = 1e-4;
end

ax = gca;
yt = get(ax,'YTick')';
dy = diff(yt);

if min(dy) > t
    return
end

if max(yt) > 0 && min(yt) < 0
    return
end

drawnow
yl = get(ax,'YTickLabels');

if strcmp(yl{1}(1),'$')
    return
end
for i = 1:length(yl)
    lens(i) = length(yl{i});
end
[~,idx] = min(lens);

dym = min(dy);
ed = floor(log10(dym));
for i = 1:length(yl)
    if i == idx
        continue
    end
    if yt(i)-yt(idx) > 0
        sign = '\hspace{-0.25em}+\hspace{-0.25em}';
        dy = floor((yt(i)-yt(idx))/dym);
    else
        sign = '\hspace{-0.25em}-\hspace{-0.25em}';
        dy = floor((yt(idx)-yt(i))/dym);  
    end
    %fprintf('Original: %s\n',yl{i});
    if dy == 1
        yl{i} = sprintf('$%s%s10^{%d}$',yl{idx},sign,ed);
    else
        yl{i} = sprintf('$%s%s%.0f\\hspace{-0.25em}\\cdot\\hspace{-0.25em} 10^{%d}$',yl{idx},sign,dy,ed);
    end
    %fprintf('New: %s\n',yl{i});
end

set(gca,'YTickLabels',yl);
drawnow

