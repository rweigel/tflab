function dock(state)

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