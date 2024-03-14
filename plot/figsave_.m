function figsave_(popts,ext)

if popts.print == 0
    return
end

if nargin == 1
    printName = popts.printOptions.printName;
else
    % Remove TeX
    ext = regexprep(ext,'\{|\}','');
    ext = regexprep(ext,'\$','');
    printName = sprintf('%s-%s',popts.printOptions.printName,ext);
end

% Need to undock window for export_fig() size commands to work.
windowstyle = get(gcf,'WindowStyle');
cf = gcf;
if ~strcmp(windowstyle,'normal')
    set(cf,'WindowStyle','normal');
    if isfield(popts,'Position')
        % Without this, setting of position does not work.
        % TODO: Find event that occurs when undock is finished and use it.
        pause(0.2);
    end
end

% Needs to be undocked for this to work.
if isfield(popts,'Position')
    set(gcf,'Position',popts.Position);
end

fname = fullfile(popts.printOptions.printDir,printName);
figsave(fname, popts.printOptions.export_fig, popts.printOptions.printFormats);

% Reset WindowStyle to initial state.
if ~strcmp(windowstyle,'normal')
    % Sometimes causes
    % Exception in thread "AWT-EventQueue-0": java.lang.NullPointerException
    set(cf,'WindowStyle','docked');
end
