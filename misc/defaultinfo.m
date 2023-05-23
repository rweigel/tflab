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

if ~iscell(S.Options.info.instr)
    if size(S.In,2) > 1
        for j = 1:size(S.In,2)
            instr{i} = sprintf('%s(:,j)',j,S.Options.info.instr);
        end
        S.Options.info.instr = instr;
    end    
end
if ~iscell(S.Options.info.outstr)
    if size(S.Out,2) > 1
        for j = 1:size(S.Out,2)
            outstr{i} = sprintf('%s(:,j)',j,S.Options.info.outstr);
        end
        S.Options.info.outstr = outstr;
    end    
end