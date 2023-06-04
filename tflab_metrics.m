function S = tflab_metrics(S)

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

In  = S.In;
Out = S.Out;

% Keep any row of Z without a NaN
% TODO: This should be handled by zinterp.
Ik = any(~isnan(S.Z),2);

% TODO: Pass zinterp options
[Zi,~] = zinterp(S.fe(Ik,:),S.Z(Ik,:),size(In,1));

OutPredicted = nan(size(Out));
for k = 1:size(In,3) % Third dimension is segment

    OutPredicted(:,:,k) = zpredict(Zi,In(:,:,k));

    S.Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    S.Metrics.MSE(1,:,k) = mse_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    S.Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k), OutPredicted(:,:,k));
 
    [S.Metrics.Coherence(:,:,k), S.Metrics.f] = ...
        coherence(Out(:,:,k), OutPredicted(:,:,k), opts, 1);

    S.Metrics.SN(:,:,k) = ...
        signaltoerror(Out(:,:,k), OutPredicted(:,:,k)-Out(:,:,k), opts, 1);
        
end

S.Out_.Predicted = OutPredicted;
