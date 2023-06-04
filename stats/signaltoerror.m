function [se,f] = signaltoerror(sig, err, opts, smoothed)

if nargin < 3
    [dftsig, f] = fftu(sig);
    dfterr = fftu(err);
    se = abs(dftsig).^2./abs(dfterr).^2;
    return
end

[dftsig,f] = dftsegments(sig, opts);
dfterr = dftsegments(err, opts);

if ~exist('smoothed', 'var')
    smoothed = 0;
end
if smoothed == 1
    w = dftweights(f, dftsig, dfterr, opts);
end

for s = 1:length(f)
    if smoothed  == 0
        se{s,1} = abs(dftsig{s}).^2./abs(dfterr{s}).^2;
    else
        se(s,:) = abs(sum(w{s}.*dftsig{s})).^2./abs(sum(w{s}.*dfterr{s})).^2;
    end
end
