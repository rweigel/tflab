function figsave(filename, opts, fmt)
%FIGSAVE Wrapper to exportgraphics
%
%   Sets background color of figure to white prior to calling exportgraphics.
%
%   FIGSAVE(filename) Calls exportgraphics(filename).
%
%   FIGSAVE(filename, opts) calls exportgraphics(filename, opts{:})
%
%   Directories are created if they do not exist.
%
%   See also PRINT, exportgraphics.


if nargin > 2
    for i = 1:length(fmt)
        fname = sprintf('%s.%s',filename, fmt{i});
        figsave(fname, opts);
    end
    return
end

if endsWith(filename, '.pdf')
    %opts = {opts{:},'-r','300'};
    opts = {opts{:},'ContentType','vector'};
end

fpath = fileparts(filename);
if ~isempty(fpath) && ~exist(fpath,'dir')
    mkdir(fpath);
    logmsg(sprintf('Created directory %s\n',fpath));
end

logmsg('Writing %s\n',filename);
if (strcmp(filename(end-3:end),'.pdf'))
    exportgraphics(gcf, filename, opts{:});
else
    %warning off export_fig:exportgraphics
    %export_fig(filename, opts{:});
    print('-dpng', '-r600', filename);
end
%exportgraphics(gcf, filename, opts{:});
logmsg('Wrote %s\n',filename);
