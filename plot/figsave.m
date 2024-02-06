function figsave(filename, opts, fmt)
%FIGSAVE Wrapper to export_fig
%
%   Sets background color of figure to white prior to calling export_fig.
%
%   FIGSAVE(filename) Calls EXPORT_FIG(filename).
%
%   FIGSAVE(filename, opts) calls EXPORT_FIG(filename, opts.export_fig{:})
%
%   Directories are created if they do not exist.
%
%   See also PRINT, EXPORT_FIG.


if nargin > 2
    for i = 1:length(fmt)
        fname = sprintf('%s.%s',filename, fmt{i});
        figsave(fname, opts);
    end
    return;
end

if endsWith(filename, '.png')
    opts = {opts{:},'-r','300'};
end

fpath = fileparts(filename);
if ~isempty(fpath) && ~exist(fpath,'dir')
    mkdir(fpath);
    logmsg(sprintf('Created directory %s\n',fpath));
end

% White background color
set(gcf,'color','w');
set(gcf,'defaultFigureColor',[1,1,1]); 

% Need to undock window for export_fig() size commands to work.
windowstyle = get(gcf,'WindowStyle');
cf = gcf;
if ~strcmp(windowstyle,'normal')
    set(cf,'WindowStyle','normal')
end

warning off export_fig:exportgraphics

logmsg(sprintf('Writing %s\n',filename));
export_fig(filename, opts{:});
logmsg(sprintf('Wrote %s\n',filename));

% Reset WindowStyle to initial state.
if ~strcmp(windowstyle,'normal')
    % Causes
    % Exception in thread "AWT-EventQueue-0": java.lang.NullPointerException
    set(cf,'WindowStyle','docked')
end
