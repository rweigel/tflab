function [se,f,clu,cll] = signaltoerror(sig, err, averaged, opts)

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

clu = NaN(length(f),size(sig,2));
cll = clu;
for s = 1:length(f)
    se(s,:) = sum(dftsig{s}.*conj(dftsig{s}),1)./sum(dfterr{s}.*conj(dfterr{s}),1);

    n = size(dftsig{s},1);
    if n > 10
        % MATLAB stats toolbox is required for bootstrp function, not used here.
        Nb = 100;
        V = nan(Nb,1);
        for comp = 1:size(dftsig{s},2)
            for b = 1:Nb
                I = randsample(n,n,1);
                V(b,1) = sum(dftsig{s}(I,comp).*conj(dftsig{s}(I,comp)),1)./sum(dfterr{s}(I,1).*conj(dfterr{s}(I,1)),1);
            end
            cll(s,comp) = prctile(V,2.5);
            clu(s,comp) = prctile(V,97.5);
        end
    end
end
