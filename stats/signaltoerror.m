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

for s = 1:length(f)
    se(s,:) = sum(dftsig{s}.*conj(dftsig{s}),1)./sum(dfterr{s}.*conj(dfterr{s}),1);
end