function [se,f,cl] = signaltoerror(sig, err, averaged, opts)

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

    n = size(dftsig{s},1);
    if n > 10
        Nb = 100;
        V = nan(Nb,1);
        for b = 1:Nb
            I = randsample(n,round(0.63*n),1);
            V(b,1) = sum(dftsig{s}(I,1).*conj(dftsig{s}(I,1)),1)./sum(dfterr{s}(I,1).*conj(dfterr{s}(I,1)),1);
        end
        nl = round((1-0.68)*Nb);
        nh = round(0.68*Nb);
        V = sort(V,1);  % Sort each column
        l = V(nl,:);    % Select the nth lowest value
        u = V(nh,:);    % Select the N-n th highest value
        cl(s,:) = [l',u'];
    end
    
end
