function S = combineStructs(S1,S2,dim,depth)
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
        S = combineStructs(S,S1{s},dim,0);
    end
    return
end
if ~exist('depth','var')
    depth = 0;
end

if depth == 0 && debug
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
    if isstruct(S1.(fns{i}))
        %logmsg(sprintf('Combining elements in structure %s\n',fns{i}));
        S.(fns{i}) = combineStructs(S1.(fns{i}),S2.(fns{i}),dim,depth+1);
    else
        % Don't combine input time series.
        nocombine = dim == 2 && depth == 0;
        nocombine = nocombine && strcmp(fns{i},'In');
        nocombine = nocombine || strcmp(fns{i},'fe');
        nocombine = nocombine || strcmp(fns{i},'f');
        if nocombine == 1
            if debug
                logmsg('Not combining field %s\n',fns{i});
            end
            S.(fns{i}) = S1.(fns{i});
        else
            %logmsg(sprintf('Combining field %s\n',fns{i}));
            if ~iscell(S1.(fns{i}))
                if debug
                    msg = 'Combining matrices in field %s across dimension %d\n';
                    logmsg(msg,fns{i},dim);
                end
                S.(fns{i}) = cat(dim,S1.(fns{i}),S2.(fns{i}));
            else
                if debug
                    msg = 'Combining matrices in cell arrays in field %s across dimension %d\n';
                    logmsg(msg,fns{i},dim);
                    S1.(fns{i})
                    S2.(fns{i})
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
