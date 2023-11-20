function figsave(filename, varargin)
%FIGSAVE Wrapper to export_fig
%
%   Sets background color of figure to white prior to calling export_fig.
%
%   FIGSAVE(filename) Calls EXPORT_FIG(filename) and creates
%   directories if needed.
%
%   FIGSAVE(filename, ...) calls EXPORT_FIG(filename, ...) and
%   creates directories if needed.
%
%   See also PRINT, EXPORT_FIG.

if endsWith(filename, '.png')
    varargin = {varargin{:},'-r','300'};
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

logmsg(sprintf('Writing %s\n',filename));
export_fig(filename, varargin{:});
logmsg(sprintf('Wrote %s\n',filename));

% Reset WindowStyle to initial state.
if ~strcmp(windowstyle,'normal')
    % Causes
    % Exception in thread "AWT-EventQueue-0": java.lang.NullPointerException
    set(cf,'WindowStyle','docked')
end
