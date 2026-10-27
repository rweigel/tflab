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
    % Confidence interval half-length
    dfteb(s,:) = dftbootstrap(bands{s});
    dfte(s,:) = t*(sqrt(vr/nr) + 1j*sqrt(vi/nr));
    % Compute percent difference for error bars
    logmsg('band = %d; n = %d\n', s, nr)
    for c = 1:size(bands{s},2)
        logmsg('  component = %d\n', c)
        if ~isreal(bands{s}(:,c))
            deltar(s,c) = abs((real(dfteb(s,c)) - real(dfte(s,c))) ./ real(dfta(s,c)) * 100);
            logmsg('    real(avg) = %.4f; real(error) = %.4f; real(error_boot) = %.4f; |real(error_diff)| = %.4f\n', ...
                real(dfta(s,c)), real(dfte(s,c)), real(dfteb(s,c)), deltar(s,c))
            deltai(s,c) = abs((imag(dfteb(s,c)) - imag(dfte(s,c))) ./ imag(dfta(s,c)) * 100);
            logmsg('    imag(avg) = %.4f; imag(error) = %.4f; imag(error_boot) = %.4f; imag(error_diff) = %.4f\n', ...
                imag(dfta(s,c)), imag(dfte(s,c)), imag(dfteb(s,c)), deltai(s,c))
        else
            delta(s,c) = abs((abs(dfteb(s,c)) - abs(dfte(s,c))) ./ abs(dfta(s,c)) * 100);
            logmsg('    avg = %.4f; error = %.4f; error_boot = %.4f; error_diff = %.4f; |error_diff| = %.4f\n', dfta(s,c), dfte(s,c), dfteb(s,c), delta(s,c), abs(dfteb(s,c)) - abs(dfte(s,c)))
        end
    end
end


if exist('deltar','var')
    logmsg('Max percent difference across all bands: real = %.4f; imag = %.4f; abs = %.4f\n', max(deltar(:)), max(deltai(:)), max(delta(:)))
else
    logmsg('Max percent difference across all bands: abs = %.4f\n', max(delta(:)))
end

end


function [delta] = dftbootstrap(band)
    %DFTBOOTSTRAP Compute bootstrap confidence intervals for DFT segments
    delta = nan(1, size(band,2));
    if size(band,1) < 20
        return
    end
    for c = 1:size(band,2) % Columns are components
        lims(c,:) = bootcl(band(:,c), @mean, 100, [2.5 97.5]);
        % Ideally we would keep upper and lower vals
        % Put warning if too much asymmetry?
        delta(c) = (lims(c,2)-lims(c,1))/2;
    end
end

