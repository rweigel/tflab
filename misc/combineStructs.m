function S = combineStructs(S1,S2,dim,ancestors)
%combineStructs Combine TFlab structures or segments

debug = 0;

if nargin == 2
    S = S1{1};
    dim = S2;
    if debug
        msg = 'Combining relevent matrices in TFLab segments across dimension %d\n';
        logmsg(msg,dim);
    end
    for s = 2:length(S1)
        % Combine segments
        S = combineStructs(S,S1{s},dim);
    end
    return
end

if ~exist('ancestors','var')
    ancestors = {};
end

if length(ancestors) == 0 && debug
    logmsg('Combining TFLab structures across dimension %d\n',dim);
end

S = struct();
fns = fieldnames(S1);
for i = 1:length(fns)
    if strcmp(fns{i},'Options')
        % Does not contain anything that can be concatenated and
        % if field exists, is same for S1 and S2.
        continue
    end
    ancestors_str = '';
    if ~isempty(ancestors)
        tmp = join(ancestors,'.');
        ancestors_str = sprintf('%s.',tmp{1});
    end
    if isstruct(S1.(fns{i}))
        ancestors{end+1} = fns{i};
        if strcmp(fns{i},'In_')
            if debug
                logmsg('Not combining elements in structure %s%s\n',ancestors_str,fns{i});
            end
            S.(fns{i}) = S1.(fns{i});
        else
            if debug
                logmsg('Combining elements in structure %s%s\n',ancestors_str,fns{i});
            end
            S.(fns{i}) = combineStructs(S1.(fns{i}),S2.(fns{i}),dim,ancestors);
        end
    else
        % Don't combine input time series.
        combine = 1;
        if (dim == 2)
            % Don't combine under certain conditions
            if strcmp(fns{i},'In')
                if isempty(ancestors)
                    % Top-level In
                    combine = 0;
                else
                    if strcmp(ancestors{end},'DFT')
                        combine = 0;
                    end
                end
            end
        end
        if strcmp(fns{i},'fe') || strcmp(fns{i},'f')
            combine = 0;
        end
        if combine == 0
            if debug
                logmsg('Not combining field %s%s\n',ancestors_str,fns{i});
            end
            S.(fns{i}) = S1.(fns{i});
        else
            if debug
                logmsg('Combining field %s%s\n',ancestors_str,fns{i});
            end
            if ~iscell(S1.(fns{i}))
                if debug
                    msg = ' Combining matrices in field %s%s across dimension %d\n';
                    logmsg(msg,ancestors_str,fns{i},dim);
                end
                S.(fns{i}) = cat(dim,S1.(fns{i}),S2.(fns{i}));
            else
                if debug
                    msg = ' Combining matrices in cell arrays in field %s%s across dimension %d\n';
                    logmsg(msg,ancestors_str,fns{i},dim);
                    %S1.(fns{i})
                    %S2.(fns{i})
                end
                for r = 1:size(S1.(fns{i}),1)
                    S.(fns{i}){r,1} = cat(dim,S1.(fns{i}){r,1},S2.(fns{i}){r,1});
                end
            end
        end
    end
end

% Use first structure's Options field. It will be same for S1 and S2.
% If Options field does not exist, not at top-level.
if isfield(S1,'Options')
    S.Options = S1.Options;
end
