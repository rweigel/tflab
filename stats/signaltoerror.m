function [se,f,clu,cll] = signaltoerror(sig, err, averaged, opts)

if ~exist('averaged', 'var')
    averaged = 0;
end

if nargin < 3 || averaged == 0
    [dftsig, f] = fftu(sig);
    dfterr = fftu(err);
    se = abs(dftsig).^2./abs(dfterr).^2;
    return
end

[dftsig,f] = dftbands(sig, opts);
dfterr = dftbands(err, opts);

se = NaN(length(f),size(sig,2));
clu = NaN(length(f),size(sig,2));
cll = clu;

p = normcdf(-2);
p = 100*[p, 1-p];

for s = 1:length(f)
    se(s,1) = sig2err(dftsig{s}, dfterr{s});
    n = size(dftsig{s},1);
    if n > 10
        % MATLAB stats toolbox is required for bootstrp function, not used here.
        Nb = 100;
        V = nan(Nb,1);
        for b = 1:Nb
            I = randsample(n,n,1);
            V(b,1) = sig2err(dftsig{s}(I), dfterr{s}(I));
        end
        cll(s,1) = prctile(V,p(1));
        clu(s,1) = prctile(V,p(2));
    end
end
end

function se = sig2err(dftsig,dfterr)
    se = sum(dftsig.*conj(dftsig))/sum(dfterr.*conj(dfterr));
end