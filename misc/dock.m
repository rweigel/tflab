function dock(state)
%DOCK Dock all figure windows or toggle dock style
%
%    DOCK() or DOCK('all') docks all windows.
%
%    DOCK('none') un-dock all windows.
%
%    DOCK('on') sets the default state for all new figures to 'docked'.
%
%    DOCK('off') sets the default state for all new figures to 'normal'.
%

if nargin == 0
    state = 'all';
end

states = {'all','none','on','off'};
if ~any(strcmp(states,state))
    tmp = sprintf('%s, ',states{:});
    error('state must be one of %s\n',tmp(1:end-2));
end

if strcmp(state,'on')
    set(0,'DefaultFigureWindowStyle','docked');
end

if strcmp(state,'off')
    set(0,'DefaultFigureWindowStyle','normal');
end

if strcmp(state,'all')
    h = findall(groot,'Type','figure');
    for f = 1:length(h)
        set(f,'WindowStyle','docked')
    end
end

if strcmp(state,'none')
    h = findall(groot,'Type','figure');
    for f = 1:length(h)
        set(f,'WindowStyle','normal')
    end
end