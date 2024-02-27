function tf = isfield_(S,fnames)
%ISFIELD_  ISFIELD but field name may contain subfields
%
%   If S.a.b = 1,
%
%   ISFIELD_(S,'a') returns 1
%   ISFIELD_(S,'a.b') returns 1
%
%   ISFIELD_(S,'x') returns 0
%   ISFIELD_(S,'x.a') returns 0
%   ISFIELD_(S,'a.x') returns 0
%   ISFIELD_(S,'a.b.x') returns 0

tf = 0;
fnamesc = split(fnames,'.');
if length(fnamesc) == 1
    tf = isfield(S,fnamesc{1});
    return
end

newfname = join(cat(1,fnamesc(2:end)),'.');
if isfield(S,fnamesc{1})
    if isfield_(S.(fnamesc{1}),newfname{1})
        tf = 1;
        return
    end
end

