function S = transferfnFD_metrics(S,opts,I,update)

if nargin < 4
    update = 0;
end

if nargin > 2 && ~isempty(I)
    % Create segments.
    if isfield(S,'Segment')
        if opts.transferfnFD.loglevel > 0
            logmsg('Segment field already exists. Will replace.\n');
        end
        % TODO: Should check if S.Segment.Intervals matches I and
        % only re-do calculation if they are not equal.s
        S = rmfield(S,'Segment');
    end
    S.Segment = struct();
    % The I variable allows start and end of segments to be removed.
    % The following assumes I(2,1,k)-I(1,1,k) is the same for all k.
    for k = 1:size(I,3)
        S.Segment.In(:,:,k) = S.In(I(1,1,k):I(2,1,k),:);
        S.Segment.Out(:,:,k) = S.Out(I(1,1,k):I(2,1,k),:);        
    end
end

if update == 1
    if isfield(S,'Metrics')
        S = rmfield(S,'Metrics');
    end
    if isfield(S,'Segment') && isfield(S.Segment,'Metrics')
        S.Segment = rmfield(S.Segment,'Metrics');
    end
end

if isfield(S,'Metrics') && isfield(S,'Segment') && isfield(S.Segment,'Metrics')
    logmsg('No metrics calculated because full and segment metrics already computed.\n');
    return
end

if ~isfield(S,'Metrics')
    % If no S.Metrics, compute
    if opts.transferfnFD.loglevel
        logmsg('Computing metrics for single segment.\n');
    end
    In = S.In;
    Out = S.Out;
else
    % If no S.Segment.Metrics, compute
    if isfield(S,'Segment') && ~isfield(S.Segment,'Metrics')
        if opts.transferfnFD.loglevel
            logmsg('Computing metrics for %d segments.\n',size(S.Segment.Out,3));
        end
        In = S.Segment.In;
        Out = S.Segment.Out;
    end
end

N = size(In,1);
Metrics = struct();

if isfield(S,'DFT') && length(S.fe) < length(S.DFT.f{1})
    Metrics.PSD = struct('Raw',struct(),'Smoothed',struct());
else
    Metrics.PSD = struct('Raw',struct());
end
Metrics.Predicted = [];

% Keep any row of Z without a NaN
% TODO: Report number of NaNs removed.
Ik = any(~isnan(S.Z),2);

% TODO: Pass zinterp options
[Zi,~] = zinterp(S.fe(Ik,:),S.Z(Ik,:),size(In,1));

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

for k = 1:size(Out,3) % Loop over segments

    Metrics.Predicted(:,:,k) = zpredict(Zi,In(:,:,k));

    Metrics.PE(1,:,k)  = pe_nonflag(Out(:,:,k),Metrics.Predicted(:,:,k));
    Metrics.MSE(1,:,k) = mse(Out(:,:,k),Metrics.Predicted(:,:,k));
    Metrics.CC(1,:,k)  = cc_nonflag(Out(:,:,k),Metrics.Predicted(:,:,k));
    
    [Metrics.PSD.Raw.In(:,:,k),Metrics.DFT.Raw.In(:,:,k),fe] = psd(In(:,:,k));
    [Metrics.PSD.Raw.Out(:,:,k),Metrics.DFT.Raw.Out(:,:,k)] = psd(Out(:,:,k));
    [Metrics.PSD.Raw.Error(:,:,k),Metrics.DFT.Raw.Error(:,:,k)] = psd(Out(:,:,k) - Metrics.Predicted(:,:,k));
    [Metrics.PSD.Raw.Predicted(:,:,k),Metrics.DFT.Raw.Predicted(:,:,k)] = psd(Metrics.Predicted(:,:,k));

    Metrics.PSD.Raw.fe = fe;
    Metrics.DFT.Raw.fe = fe;
    
    Metrics.SN.Raw(:,:,k)  = Metrics.PSD.Raw.Out(:,:,k)./Metrics.PSD.Raw.Error(:,:,k);
    Metrics.Coherence.Raw(:,:,k) = coherence(Out(:,:,k),Metrics.Predicted(:,:,k));

    if smoothed
        [Metrics.PSD.Smoothed.In(:,:,k),Metrics.DFT.Smoothed.In(:,:,k),fe] = psd(In(:,:,k),opts,N);
        [Metrics.PSD.Smoothed.Out(:,:,k),Metrics.DFT.Smoothed.Out(:,:,k)] = psd(Out(:,:,k),opts,N);
        [Metrics.PSD.Smoothed.Error(:,:,k),Metrics.DFT.Smoothed.Error(:,:,k)] = psd(Out(:,:,k) - Metrics.Predicted(:,:,k),opts,N);
        [Metrics.PSD.Smoothed.Predicted(:,:,k),Metrics.DFT.Smoothed.Predicted(:,:,k)] = psd(Metrics.Predicted(:,:,k),opts,N);

        Metrics.PSD.Smoothed.fe = fe;
        Metrics.DFT.Smoothed.fe = fe;

        Metrics.SN.Smoothed(:,:,k)  = Metrics.PSD.Smoothed.Out(:,:,k)./Metrics.PSD.Smoothed.Error(:,:,k);
        Metrics.Coherence.Smoothed(:,:,k) = coherence(Out(:,:,k),Metrics.Predicted(:,:,k),opts,N);
    end

    filts = {'Window','Prewhiten','Zeropad'};
    for f=1:length(filts) 
        filt = filts{f};
        if isfield(S,filt)
            S.(filt).PSD = struct();
                [S.(filt).PSD.Raw.In,S.(filt).DFT.Raw.In,fe] = psd(S.(filt).In);
                [S.(filt).PSD.Raw.Out,S.(filt).DFT.Raw.Out] = psd(S.(filt).Out);
                S.(filt).PSD.Raw.fe = fe;
                S.(filt).DFT.Raw.fe = fe;
                [S.(filt).PSD.Smoothed.In,S.(filt).DFT.In,fe] = psd(S.(filt).In,opts,N);
                [S.(filt).PSD.Smoothed.Out,S.(filt).DFT.Smoothed.Out] = psd(S.(filt).Out,opts,N);
                S.(filt).PSD.Smoothed.fe = fe;
                S.(filt).DFT.Smoothed.fe = fe;
        end
    end
    
end

if ~isfield(S,'Metrics')
    S.Metrics = Metrics;
else
    S.Segment.Metrics = Metrics;
end

if isfield(S,'Segment') && ~isfield(S.Segment,'Metrics')
    % If no S.Metrics and no S.Segment.Metrics, then 
    % S.Metrics was computed first. This catches case where both metrics
    % need to be calculated.
    S.Segment = transferfnFD_metrics(S,opts);
    return
end