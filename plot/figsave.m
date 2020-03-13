function figsave(y,fname)

if ~y
    return
end

addpath([fileparts(mfilename('fullpath')),'/export_fig']);

[ST, I] = dbstack('-completenames', 1);

if length(ST) > 0
    fprintf('%s.m: Writing %s\n',ST(1).name,fname);
else
    fprintf('Writing %s\n',fname);
end

set(gcf,'color','w');
%export_fig(fname,'-transparent');
export_fig(fname);

if length(ST) > 0
    fprintf('%s.m: Wrote %s\n',ST(1).name,fname);
else
    fprintf('Wrote %s\n',fname);
end
    
