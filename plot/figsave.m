function figsave(popts,default)

addpath([fileparts(mfilename('fullpath')),'/../misc']); % logmsg.m
addpath([fileparts(mfilename('fullpath')),'/export_fig']);

if isfield(popts,'filename') && ~isempty(popts.filename)
    filename = popts.filename;
else
    if nargin < 2
        % Use calling script filename.
        filename = caller();
        if isempty(filename)
            % Will happen if figsave() called from command line.
            filename = 'figsave';
        end
    else
        filename = default;
    end
end

fpath = fileparts(filename);
if ~isempty(fpath) && ~exist(fpath,'dir')
    mkdir(fpath);
    logmsg(sprintf('Created directory %s\n',fpath));
end

set(gcf,'color','w');
set(gcf,'defaultFigureColor', [1,1,1]); % Background color to white.

fns = fieldnames(popts.savefmt);
for i = 1:length(fns)
    if popts.savefmt.(fns{i})
        fname = [filename,'.',fns{i}];
        logmsg(sprintf('Writing %s\n',fname));
        export_fig(fname);
        logmsg(sprintf('Wrote %s\n',fname));
    end
end