function dfta = dftaverage(band, weights)
%DFTAVERAGE

for s = 1:length(band)
    nc = size(band{s},2);
    nr = size(band{s},1);
    if ~exist('weights','var')
        dfta(s,:) = sum(band{s},1)/nr;
    else
        dfta(s,:) = sum(weights{s}.*band{s},1);
    end
end


