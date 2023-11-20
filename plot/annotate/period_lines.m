function period_lines(show_labels)
%PERIOD_LINES Vertical lines and labels at periods in current x-axis range.
%
%   Usage PERIOD_LINES()

% TODO: Allow timeunit to be passed.

if nargin < 1
    show_labels = 1;
end

drawnow; % Updates axes limits

held = ishold();
if ~held
    hold on;
end

xl = get(gca,'XLim');
yl = get(gca,'YLim');

at = [60,   60*60, 60*60*6, 60*60*12, 60*60*24, 60*60*24*5, 60*60*24*10, 60*60*24*30];
ls = {'1m','1h' , '6h',    '12h',    '1d',     '5d',       '10d',       '30d'};
%at = at*freqsf;

for i = 1:length(at)
    if at(i) >= xl(end)
        break
    end
    plot([at(i),at(i)],yl,'--','Color',[0.5,0.5,0.5]);
    if show_labels == 1
        text(at(i),yl(1),ls{i},'VerticalAlignment','bottom');
    end
end

if ~held
    % Reset hold state to initial.
    hold off;
end
