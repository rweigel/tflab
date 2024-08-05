function [cxy,f] = coherence(x, y, averaged, opts)

if ~exist('averaged', 'var')
    averaged = 0;
end

if nargin < 3 || averaged == 0
    [cxy,f] = mscohere(x,y);
    return
end

[dftsegsx,f,fe] = dftbands(x, opts);
dftsegsy = dftbands(y, opts);
[wx,wy] = dftweights(f, dftsegsx, dftsegsy, opts);

cxy = NaN(length(f),1);
clu = NaN(length(f),1);
cll = clu;

p = normcdf(-2);
p = 100*[p, 1-p];

for s = 1:length(f)
    if length(dftsegsx{s}) == 1
        % If only one DFT value, coherence is 1 by definition. This
        % avoids cxy = NaN when either dftsegsx{s} or dftsegsy{s} are
        % identically zero.
        cxy(s,1) = 1;
        continue
    end
    cxy(s,1) = coh(dftsegsx{s}, dftsegsy{s});
    n = size(dftsegsx{s},1);
    if n > 10
        % MATLAB stats toolbox is required for bootstrp function, not used here.
        Nb = 100;
        V = nan(Nb,1);
        for b = 1:Nb
            I = randsample(n,n,1);
            V(b,1) = coh(dftsegsx{s}(I), dftsegsy{s}(I));
        end
        cll(s,1) = prctile(V,p(1));
        clu(s,1) = prctile(V,p(2));
    end
end
f = fe;
end

function cxy = coh(x, y)
    sxx = abs(sum(x.*conj(x)));
    syy = abs(sum(y.*conj(y)));
    sxy = abs(sum(conj(x).*y));
    cxy = sqrt(sxy/sqrt(sxx*syy));
end

