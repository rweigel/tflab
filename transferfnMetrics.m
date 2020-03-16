function S = transferfnMetrics(S,opts,I)
%TRANSFERFNMETRICS

if nargin == 3
    % Create segments.
    % The following assumes I(2,1,k)-I(1,1,k) is constant.
    if isfield(S,'Segment')
        logmsg('Segment field already exists. Will replace.\n');
        S = rmfield(S,'Segment');
    end
    S.Segment = struct();
    for k = 1:size(I,3)
        S.Segment.In(:,:,k) = S.In(I(1,1,k):I(2,1,k),:);
        S.Segment.Out(:,:,k) = S.Out(I(1,1,k):I(2,1,k),:);        
    end
end

if isfield(S,'Metrics') && isfield(S,'Segment') && isfield(S.Segment,'Metrics')
    logmsg('All metrics computed already. No metrics calculated.\n');
    return
end

if ~isfield(S,'Metrics')
    % If no S.Metrics, compute
    In = S.In;
    Out = S.Out;
else
    % If no S.Segment.Metrics, compute
    if isfield(S,'Segment') && ~isfield(S.Segment,'Metrics')
        In = S.Segment.In;
        Out = S.Segment.Out;
    end
end

if isfield(S,'Segment')
    % In this case, S.In will be a single segment (the full time series).
    % This forces PSD calculations on full time series to have evaluation
    % frequencies that match those for the segments.
    N = size(S.Segment.In,1);
else

end

N = size(In,1);
Metrics = struct();
Metrics.PSD = struct();
Metrics.Predicted = [];

[Zi,fi] = zinterp(S.fe,S.Z,size(In,1));

for k = 1:size(Out,3)

    Metrics.Predicted(:,:,k) = zpredict(Zi,In(:,:,k));

    [Metrics.PSD.In(:,:,k),Metrics.fe] = smoothSpectra(In(:,:,k),opts,N);
    Metrics.PSD.Out(:,:,k)   = smoothSpectra(Out(:,:,k),opts,N);

    Metrics.PSD.Error(:,:,k) = smoothSpectra(Out(:,:,k) - Metrics.Predicted(:,:,k),opts,N);
    Metrics.PSD.Predicted(:,:,k) = smoothSpectra(Metrics.Predicted(:,:,k),opts,N);

    Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k),Metrics.Predicted(:,:,k));
    Metrics.MSE(1,:,k) = mse(Out(:,:,k),Metrics.Predicted(:,:,k));
    Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k),Metrics.Predicted(:,:,k));
    Metrics.SN(:,:,k)  = Metrics.PSD.Out(:,:,k)./Metrics.PSD.Error(:,:,k);
    Metrics.Coherence(:,:,k) = smoothCoherence(...
                                    Out(:,:,k),...
                                    Metrics.Predicted(:,:,k),opts,N);
end

if ~isfield(S,'Metrics')
    S.Metrics = Metrics;
else
    S.Segment.Metrics = Metrics;
end

if isfield(S,'Segment') && ~isfield(S.Segment,'Metrics')
    % If no metrics in S and no metrics in S.Segment, then 
    % S.Metrics was computed first. This catches case where both metrics
    % need to be calculated.
    S.Segment = transferfnMetrics(S,opts);
    return
end