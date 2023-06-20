function S = combineStructs(S1,S2,dim)
%combineStructs Combine tflab structures

    if nargin == 2
        S = S1{1};
        dim = S2;
        for s = 2:length(S1)
            % Combine segments
            S = combineStructs(S,S1{s},dim);
        end
        return;
    end

    S = struct();    
    fns = fieldnames(S1);
    for i = 1:length(fns)
        if strcmp(fns{i},'Options')
            % Does not contain anything that can be concatenated and
            % if field exists, is same for S1 and S2.
            continue;
        end
        if isstruct(S1.(fns{i}))
            S.(fns{i}) = combineStructs(S1.(fns{i}),S2.(fns{i}),dim);
        else
            if (dim == 2 && strcmp(fns{i},'In')) || strcmp(fns{i},'fe')
                S.(fns{i}) = S1.(fns{i});
            else
                S.(fns{i}) = cat(dim,S1.(fns{i}),S2.(fns{i}));
            end
        end
    end
    % Use first structure's Options field. It will be same for S1 and S2.
    % If Options field does not exist, not at top-level.
    if isfield(S1,'Options')
        S.Options = S1.Options;
    end

end
