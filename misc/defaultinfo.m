function S = defaultinfo(S)

if iscell(S)
    for s = 1:length(S)
        S{s} = defaultinfo(S{s});
    end
    return;
end

dopts = tflab_options();

if ~isfield(S,'Options')
    S.Options = struct();
end
if ~isfield(S.Options,'description')
    S.Options.description = '';
end
if ~isfield(S.Options,'info')
    S.Options.info = dopts.info;
end

S.Options.info = updatedefaults(dopts.info,S.Options.info);


