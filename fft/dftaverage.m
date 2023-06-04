function dfta = dftaverage(dftsegments, weights)

for s = 1:length(dftsegments)
    nc = size(dftsegments{s},2);
    nr = size(dftsegments{s},1);
    if ~exist('weights','var')
        dfta(s,:) = sum(dftsegments{s},1)/nr;
    else
        w = repmat(weights{s},1,nc);
        dfta(s,:) = sum(w.*dftsegments{s},1);
    end
end


