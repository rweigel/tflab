function S = defaultinfo(S)

dopts = tflab_options();

if ~isfield(S,'Options')
    S.Options = struct();
end
if ~isfield(S.Options,'description')
    S.Options.description = '';
end
if ~isfield(S.Options,'info')
    S.Options.info = dopts.info;
    return;
end

S.Options.info = updatedefaults(dopts.info,S.Options.info);
