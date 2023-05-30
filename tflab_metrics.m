function S = tflab_metrics(S)

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

In  = S.In;
Out = S.Out;

% Keep any row of Z without a NaN
% TODO: Report number of NaNs removed.
Ik = any(~isnan(S.Z),2);

% TODO: Pass zinterp options
[Zi,~] = zinterp(S.fe(Ik,:),S.Z(Ik,:),size(In,1));

PSD = struct();
for k = 1:size(In,3)
    S.Metrics.Predicted(:,:,k) = zpredict(Zi,In(:,:,k));

    Predicted = S.Metrics.Predicted;

    S.Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k), Predicted(:,:,k));
    S.Metrics.MSE(1,:,k) = mse(Out(:,:,k), Predicted(:,:,k));
    S.Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k), Predicted(:,:,k));

    [S.Metrics.PSD.Raw.Error(:,:,k), S.Metrics.DFT.Raw.Error(:,:,k),fe] = psd(Out(:,:,k) - Predicted(:,:,k));
    [S.Metrics.PSD.Raw.Predicted(:,:,k), S.Metrics.DFT.Raw.Predicted(:,:,k)] = psd(Predicted(:,:,k));

    S.Metrics.PSD.Raw.fe = fe;
    S.Metrics.DFT.Raw.fe = fe;

    [S.Metrics.PSD.Raw.In(:,:,k), S.Metrics.DFT.Raw.In(:,:,k), S.Metrics.DFT.Raw.fe] = psd(In(:,:,k));
    [S.Metrics.PSD.Raw.Out(:,:,k), S.Metrics.DFT.Raw.Out(:,:,k)] = psd(Out(:,:,k));

    S.Metrics.SN.Raw(:,:,k) = S.Metrics.PSD.Raw.Out(:,:,k)./S.Metrics.PSD.Raw.Error(:,:,k);
    S.Metrics.Coherence.Raw(:,:,k) = coherence(Out(:,:,k), Predicted(:,:,k));

    fnargs = opts.fd.evalfreq.functionargs;
    smoothed = 1;
    if length(fnargs) == 3 && strcmp('linear',fnargs{2})
        if length(fnargs{1}) > 1 && fnargs{1}(2) == 0
            % evalfreq(N,[dN,0],'linear')
            smoothed = 0;
        end
        if length(fnargs{1}) == 1 && fnargs{1}(1) == 0
            % evalfreq(N,0,'linear')
            smoothed = 0;
        end
    end

    if smoothed
        [S.Metrics.PSD.Smoothed.In(:,:,k), S.Metrics.DFT.Smoothed.In(:,:,k), fe] = psd(In(:,:,k), opts);
        [S.Metrics.PSD.Smoothed.Out(:,:,k), S.Metrics.DFT.Smoothed.Out(:,:,k)] = psd(Out(:,:,k), opts);
        [S.Metrics.PSD.Smoothed.Error(:,:,k), S.Metrics.DFT.Smoothed.Error(:,:,k)] = psd(Out(:,:,k) - S.Metrics.Predicted(:,:,k), opts);
        [S.Metrics.PSD.Smoothed.Predicted(:,:,k), S.Metrics.DFT.Smoothed.Predicted(:,:,k)] = psd(S.Metrics.Predicted(:,:,k), opts);

        S.Metrics.PSD.Smoothed.fe = fe;
        S.Metrics.DFT.Smoothed.fe = fe;

        S.Metrics.SN.Smoothed(:,:,k) = S.Metrics.PSD.Smoothed.Out(:,:,k)./S.Metrics.PSD.Smoothed.Error(:,:,k);
        S.Metrics.Coherence.Smoothed(:,:,k) = coherence(Out, S.Metrics.Predicted(:,:,k),opts);
    end

    filts = {'Window','Prewhiten','Zeropad','Detrend'};
    for f = 1:length(filts) 
        filt = filts{f};
        if isfield(S,filt)
            S.(filt).PSD = struct();
                [S.(filt).PSD.Raw.In,S.(filt).DFT.Raw.In,fe] = psd(S.(filt).In);
                [S.(filt).PSD.Raw.Out,S.(filt).DFT.Raw.Out] = psd(S.(filt).Out);
                S.(filt).PSD.Raw.fe = fe;
                S.(filt).DFT.Raw.fe = fe;
                [S.(filt).PSD.Smoothed.In,S.(filt).DFT.In,fe] = psd(S.(filt).In,opts);
                [S.(filt).PSD.Smoothed.Out,S.(filt).DFT.Smoothed.Out] = psd(S.(filt).Out,opts);
                S.(filt).PSD.Smoothed.fe = fe;
                S.(filt).DFT.Smoothed.fe = fe;
        end
    end
end
