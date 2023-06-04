function [cxy,f] = coherence(x, y, opts, smoothed)

if nargin < 3
    [dftx, f] = fftu(x);
    dfty = fftu(y);
    dftxy = conj(dftx).*dfty;
    cxy = abs(dftxy)./(abs(dftx).*abs(dfty));
    return
end

[dftsegsx,f] = dftsegments(x, opts);
dftsegsy = dftsegments(y, opts);

if ~exist('smoothed', 'var')
    smoothed = 0;
end
if smoothed == 1
    w = dftweights(f, dftsegsx, dftsegsy, opts);
end

for s = 1:length(f)
    dftx = dftsegsx{s};
    dfty = dftsegsy{s};
    dftxy  = dftx.*conj(dfty);
    if smoothed == 0       
        cxy{s,1} = abs(dftxy)/(abs(dftx).*abs(dfty));
    else
        dftx  = sum(w{s}.*dftsegsx{s});
        dfty  = sum(w{s}.*dftsegsy{s});
        dftxy = sum((w{s}.*dftsegsx{s}).*conj(w{s}.*dftsegsy{s}));
        cxy(s,:) = abs(dftxy)/(abs(dftx).*abs(dfty));
    end
end
