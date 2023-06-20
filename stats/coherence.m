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
        % If only one dft value, coherence is 1 by definition. This
        % avoids cxy = NaN when either dftsegsx{s} or dftsegsy{s} are
        % identically zero.
        cxy(s,:) = 1;
        continue;
    end
    sxx = sum(abs(wx{s}.*dftsegsx{s}).^2,1);
    syy = sum(abs(wy{s}.*dftsegsy{s}).^2,1);
    sxy = sum(abs( conj(wx{s}.*dftsegsx{s}).*wy{s}.*dftsegsy{s} ),1);
    cxy(s,:) = sxy.^2/(sxx*syy);
end

f = fe;
