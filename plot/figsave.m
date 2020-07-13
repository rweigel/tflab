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

addpath([fileparts(mfilename('fullpath')),'/../misc']); % logmsg.m
addpath([fileparts(mfilename('fullpath')),'/export_fig']);

fpath = fileparts(filename);
if ~isempty(fpath) && ~exist(fpath,'dir')
    mkdir(fpath);
    logmsg(sprintf('Created directory %s\n',fpath));
end

% White background color
set(gcf,'color','w');
set(gcf,'defaultFigureColor',[1,1,1]); 

logmsg(sprintf('Writing %s\n',filename));
export_fig(filename, varargin{:});
logmsg(sprintf('Wrote %s\n',filename));
