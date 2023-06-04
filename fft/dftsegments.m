function [segs, f, fe] = dftsegments(x, opts)

[fe,Ic,Ne] = opts.fd.evalfreq.function(...
                size(x,1),opts.fd.evalfreq.functionargs{:});

[dftu, fu] = fftu(x);

for j = 1:length(Ic)
    r = Ic(j)-Ne(j):Ic(j)+Ne(j); % Index range
    segs{j,1} = dftu(r,:);
    f{j,1} = fu(r,:);
end
