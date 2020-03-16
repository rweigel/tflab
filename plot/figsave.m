function figsave(y,fname)

if ~y
    return
end

addpath([fileparts(mfilename('fullpath')),'/export_fig']);

logmsg(sprintf('Writing %s\n',fname));

set(gcf,'color','w');
set(gcf,'defaultFigureColor', [1,1,1]); % Background color to white.

% Save figure
%export_fig(fname,'-transparent');
export_fig(fname);

logmsg(sprintf('Wrote %s\n',fname));