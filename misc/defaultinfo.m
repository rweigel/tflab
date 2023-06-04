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

if size(S.Out,2) == 1
    if size(S.In,2) == 1
        zstrs = {'Z'};
        rhostrs = {'\rho^a'};
        phistrs = {'\phi'};
    end
    if size(S.In,2) == 2
        zstrs = {'Z_{x}','Z_{y}'};
        rhostrs = {'\rho^a_{x}','\rho^a_{y}'};
        phistrs = {'\phi_{x}','\phi_{y}'};
    end
    if size(S.In,2) == 3
        zstrs = {'Z_{x}','Z_{y}','Z_{z}'};
        rhostrs = {'\rho^a_{x}','\rho^a_{y}','\rho^a_{z}'};
        phistrs = {'\phi_{x}','\phi_{y}','\phi_{z}'};
    end
end
if size(S.Out,2) == 2 && size(S.In,2) == 2
    zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{yx}','\rho^a_{yy}'};
    phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
end
if exist('zstrs','var')
    % TODO: Handle other cases for dims of In and Out
    if ~isfield(S.Options.info,'zstrs')
        S.Options.info.zstrs = zstrs;
    end
    if ~isfield(S.Options.info,'rhostrs')
        S.Options.info.rhostrs = rhostrs;
    end
    if ~isfield(S.Options.info,'phistrs')
        S.Options.info.phistrs = phistrs;
    end
end

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
