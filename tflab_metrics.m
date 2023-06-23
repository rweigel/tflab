function S = tflab_metrics(S,onsegments)

if nargin < 2
    onsegments = 0;
end

opts = S.Options;
if onsegments
    In  = S.Segment.In;
    Out = S.Segment.Out;
else
    In  = S.In;
    Out = S.Out;
end

assert(size(S.Z,2) == size(In,2),'S.Z must have same number of columns as S.In');

[Zi,~] = zinterp(S.fe,S.Z,size(In,1));

OutPredicted = nan(size(Out));

for k = 1:size(In,3) % Third dimension is segment

    OutPredicted(:,:,k) = zpredict(Zi,In(:,:,k));

    Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    Metrics.MSE(1,:,k) = mse_nonflag(Out(:,:,k), OutPredicted(:,:,k));
    Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k), OutPredicted(:,:,k));
 
    [Metrics.Coherence(:,:,k), S.Metrics.fe] = ...
        coherence(Out(:,:,k), OutPredicted(:,:,k), 1, opts);

    Error(:,:,k) = OutPredicted(:,:,k)-Out(:,:,k);

    Metrics.SN(:,:,k) = ...
        signaltoerror(Out(:,:,k), Error(:,:,k), 1, opts);

    DFTError(:,:,k) = dftbands(Error(:,:,k), opts);    
    DFTPredicted(:,:,k) = dftbands(OutPredicted, opts);
    
end
if onsegments
    S.Segment.Metrics = Metrics;
    S.Segment.Out_.Error = Error;
    S.Segment.Out_.Predicted = OutPredicted;    
    S.Segment.DFT.Out_.Error = DFTError;
    S.Segment.DFT.Out_.Predicted = DFTPredicted;
else
    S.Metrics = Metrics;
    S.Out_.Error = Error;
    S.Out_.Predicted = OutPredicted;
    S.DFT.Out_.Error = DFTError;
    S.DFT.Out_.Predicted = DFTPredicted;
end