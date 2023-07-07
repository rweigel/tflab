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

assert(size(Out,2) == size(S.Z,2)/size(In,2),'size(Out,2) must equal size(Z,2)/size(In,2)');

[Zi,~] = zinterp(S.fe,S.Z,size(In,1));

OutPredicted = nan(size(Out));

for j = 1:size(Out,2) % Second dimension is component

    zcols = (1:size(In,2)) + (j-1)*size(In,2);

    for k = 1:size(In,3) % Third dimension is segment
 
        OutPredicted(:,j,k) = zpredict(Zi(:,zcols),In(:,:,k));
        Metrics.PE(1,j,k)  = pe_nonflag(Out(:,j,k), OutPredicted(:,j,k));
        Metrics.MSE(1,j,k) = mse_nonflag(Out(:,j,k), OutPredicted(:,j,k));
        Metrics.CC(1,j,k)  = cc_nonflag(Out(:,j,k), OutPredicted(:,j,k));

        [Metrics.Coherence(:,j,k), Metrics.fe] = ...
            coherence(Out(:,j,k), OutPredicted(:,j,k), 1, opts);

        Error(:,j,k) = OutPredicted(:,j,k)-Out(:,j,k);

        Metrics.SN(:,j,k) = ...
            signaltoerror(Out(:,j,k), Error(:,j,k), 1, opts);

    end
end

for k = 1:size(In,3)
    [DFTError(:,:,k),f,fe] = dftbands(Error(:,:,k), opts);
    DFTPredicted(:,:,k) = dftbands(OutPredicted(:,:,k), opts);
end

if onsegments
    S.Segment.Metrics = Metrics;
    S.Segment.Out_.Predicted = OutPredicted;    
    S.Segment.DFT.f = f;
    S.Segment.DFT.fe = fe;
    S.Segment.DFT.Out_.Error = DFTError;
    S.Segment.DFT.Out_.Predicted = DFTPredicted;
else
    S.Metrics = Metrics;
    S.Out_.Error = Error;
    S.Out_.Predicted = OutPredicted;
    S.DFT.f = f;
    S.DFT.fe = fe;
    S.DFT.Out_.Error = DFTError;
    S.DFT.Out_.Predicted = DFTPredicted;
end