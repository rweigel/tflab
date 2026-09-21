function [dfta,dfte] = dftaverage(bands, weights)
%DFTAVERAGE

for s = 1:length(bands)
    nr = size(bands{s},1); % # of records
    if nr < 11
        continue
    end
    dfte(s,:) = nan(1, size(bands{s},2));
    t = tinv(0.975,nr-1);
    if ~exist('weights','var')
        dfta(s,:) = mean(bands{s},1);
    else
        % When weights are not all the same, error bar estimate is biased.
        % For unbiased estimate derivation, see
        % https://mathoverflow.net/a/11870 (not confirmed, but seems correct)
        dfta(s,:) = sum(weights{s}.*bands{s},1)./sum(weights{s},1);
    end
    % See page 288 of Devore 8th ed.
    % https://drive.google.com/file/d/11Ggp-RNoknu7ARu95s54hvOsQMv0AgR-/
    vr = var(real(bands{s}),0,1);
    vi = var(imag(bands{s}),0,1);
    dfte(s,:) = t*(sqrt(vr/nr) + 1j*sqrt(vi/nr));
end