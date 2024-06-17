function [cxy,f] = coherence(x, y, averaged, opts)

if ~exist('averaged', 'var')
    averaged = 0;
end

if nargin < 3 || averaged == 0
    [cxy,f] = mscohere(x,y);
    return;
end

[dftsegsx,f,fe] = dftbands(x, opts);
dftsegsy = dftbands(y, opts);
[wx,wy] = dftweights(f, dftsegsx, dftsegsy, opts);

for s = 1:length(f)
    if length(dftsegsx{s}) == 1
        % If only one DFT value, coherence is 1 by definition. This
        % avoids cxy = NaN when either dftsegsx{s} or dftsegsy{s} are
        % identically zero.
        cxy(s,:) = 1;
        continue
    end
    sxx = abs(sum(dftsegsx{s}.*conj(dftsegsx{s}),1));
    syy = abs(sum(dftsegsy{s}.*conj(dftsegsy{s}),1));
    sxy = abs(sum(conj(dftsegsx{s}).*dftsegsy{s},1));
    cxy(s,:) = sqrt(sxy/sqrt(sxx*syy));
end

f = fe;
