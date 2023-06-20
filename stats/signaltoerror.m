function [se,f] = signaltoerror(sig, err, averaged, opts)

if ~exist('averaged', 'var')
    averaged = 0;
end

if nargin < 3 || averaged == 0
    [dftsig, f] = fftu(sig);
    dfterr = fftu(err);
    se = abs(dftsig).^2./abs(dfterr).^2;
    return;
end

[dftsig,f] = dftbands(sig, opts);
dfterr = dftbands(err, opts);

w = dftweights(f, dftsig, dftsig, opts);
for s = 1:length(f)
    se(s,:) = abs(sum(w{s}.*dftsig{s},1)).^2./abs(sum(w{s}.*dfterr{s},1)).^2;
end
