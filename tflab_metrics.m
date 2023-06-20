function S = tflab_metrics(S)

opts = S.Options;

In  = S.In;
Out = S.Out;

assert(size(S.Z,2) == size(In,2),'S.Z must have same number of columns as S.In');

[Zi,~] = zinterp(S.fe,S.Z,size(In,1));

OutPredicted = nan(size(Out));

for k = 1:size(In,3) % Third dimension is segment

    OutPredicted(:,:,k) = zpredict(Zi,In(:,:,k));

    S.Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    S.Metrics.MSE(1,:,k) = mse_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    S.Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k), OutPredicted(:,:,k));
 
    [S.Metrics.Coherence(:,:,k), S.Metrics.fe] = ...
        coherence(Out(:,:,k), OutPredicted(:,:,k), 1, opts);

    Error(:,:,k) = OutPredicted(:,:,k)-Out(:,:,k);

    S.Metrics.SN(:,:,k) = ...
        signaltoerror(Out(:,:,k), Error(:,:,k), 1, opts);

    S.DFT.Out_.Error(:,:,k) = dftbands(Error(:,:,k), opts);    
    S.DFT.Out_.Predicted(:,:,k) = dftbands(OutPredicted, opts);
    
end

S.Out_.Error = Error;
S.Out_.Predicted = OutPredicted;
