function [dfta,dfte] = dftaverage(bands, weights)
%DFTAVERAGE

for s = 1:length(bands)
    nr = size(bands{s},1);
    if ~exist('weights','var')
        dfta(s,:) = mean(bands{s},1);
        dfte(s,:) = nan(1, size(bands{s},2));
        if length(bands) > 5
            % See page 288 of Devore 8th ed.
            % https://drive.google.com/file/d/11Ggp-RNoknu7ARu95s54hvOsQMv0AgR-/%E2%98%85%E2%98%85%E2%98%85%E2%98%85remove%E2%98%85%E2%98%85%E2%98%85%E2%98%85
            t = tinv(0.975,size(bands{s},1)-1);
            dfte(s,:) = t*sqrt(var(bands{s},0,1)/nr);
        end
    else
        dfta(s,:) = sum(weights{s}.*bands{s},1);
        dfte(s,:) = nan(1, size(bands{s},2));
        if length(bands) > 5
            % This is a biased estimate. For unbiased estimate derivation, see
            % https://mathoverflow.net/a/11870 (not confirmed, but seems correct)
            t = tinv(0.975,size(bands{s},1)-1);
            dfte(s,:) = t*sqrt(var(bands{s},0,1)/nr).*sum(weights{s}.^2,1);
        end
    end
end