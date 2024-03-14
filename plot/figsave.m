function figsave(filename, opts, fmt)
%FIGSAVE Wrapper to export_fig
%
%   Sets background color of figure to white prior to calling export_fig.
%
%   FIGSAVE(filename) Calls EXPORT_FIG(filename).
%
%   FIGSAVE(filename, opts) calls EXPORT_FIG(filename, opts{:})
%
%   Directories are created if they do not exist.
%
%   See also PRINT, EXPORT_FIG.


if nargin > 2
    for i = 1:length(fmt)
        fname = sprintf('%s.%s',filename, fmt{i});
        figsave(fname, opts);
    end
    return
end

if endsWith(filename, '.png')
    opts = {opts{:},'-r','300'};
end

fpath = fileparts(filename);
if ~isempty(fpath) && ~exist(fpath,'dir')
    mkdir(fpath);
    logmsg(sprintf('Created directory %s\n',fpath));
end

warning off export_fig:exportgraphics

logmsg('Writing %s\n',filename);
export_fig(filename, opts{:});
logmsg('Wrote %s\n',filename);
