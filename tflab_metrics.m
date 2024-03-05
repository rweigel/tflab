function S = tflab_metrics(S,onsegments)
%TFLAB_METRICS Compute metrics using TFLab structure
%
%   S = TFLAB_METRICS(S) adds metrics calculations to TFlab structure.
%
%   If S.In is a matrix, structures S.Metrics, S.Out_, and S.DFT.Out_
%   are computed using calculations based on S.In, S.Out, S.DFT, S.Z,
%   and S.fe.
%
%   If S has a Segment field, S.Segment.Metrics, S.Segment.Out_, and
%   S.Segment.DFT.Out_ are computed using matrices S.Segment.In,
%   S.Segment.Out, S.Segment.DFT, S.Segment.Z, and S.Segment.fe. If
%   S.Segment.Z does not exist, S.Z and S.fe are used.
%
%   If S.In is a cell array, cell array elements S.Metrics{i},
%   S.Out_{i}, and S.DFT.Out_{i} are computed for each i of S using
%   calculations based on matrices S.In{i}, S.Out{i}, S.DFT{i}, S.Z,
%   and S.fe. If S has a Segment field, metrics are computed on the 
%   segments in the same way as described in the previous paragraph.

if nargin < 2
    onsegments = 0;
end

opts = S.Options;
if onsegments == 0
    if isfield_(S,'Regression.ErrorEstimates')
        % Compute Metrics.ErrorEstimates based on content of
        % S.Regression.ErrorEstimates.Parametric
        logmsg('Computing error estimates.\n');
        Metrics.ErrorEstimates = error_estimates(S,0);
    end
    if ~iscell(S.In)
        logmsg('Computing metrics on full unsegmented data using top-level Z.\n');
        In  = S.In;
        Out = S.Out;
        DFT = S.DFT;
        Z = S.Z;
        fe = S.fe;
        dZ = [];
        if isfield(S,'dZ')
            dZ = S.dZ;
        end
    else
        % S.In is a cell array of intervals.
        Ni = length(S.In);
        msg = sprintf('Computing metrics on %d intervals using top-level Z.\n',Ni);
        logmsg(msg);
        In_intervals = S.In;
        Out_intervals = S.Out;
        DFT_intervals = S.DFT;
        for i = 1:Ni
            logmsg('Computing metrics on interval %d.\n',i);
            S.In = In_intervals{i};
            S.Out = Out_intervals{i};
            S.DFT = DFT_intervals{i};
            Ss = tflab_metrics(S,0);
            Metrics{i} = Ss.Metrics;
        end
        S.In = In_intervals;
        S.Out = Out_intervals;
        S.DFT = DFT_intervals;
        if isfield(S,'Segment')
            logmsg('Computing metrics on segments.\n')
            S = tflab_metrics(S,1);
        end
        S.Metrics = Metrics;
        return
    end
else
    assert(isfield(S,'Segment'), 'S must have a Segment field')
    In  = S.Segment.In;
    Out = S.Segment.Out;
    DFT = S.Segment.DFT;
    dZ = [];
    if isfield(S.Segment,'Z')
        logmsg('Computing metrics on segments using segment Zs.\n');
        logmsg('Computing error estimates.\n');
        S.Metrics.ErrorEstimates = error_estimates(S,1);
        Z = S.Segment.Z;
        fe = S.Segment.fe;
        if isfield(S.Segment,'dZ')
            dZ = S.Segment.dZ;
        end
    else
        % Use top-level Z to compute metrics on each segement.
        logmsg('Computing metrics on segments using top-level Z.\n');
        Z = S.Z;
        fe = S.fe;
        if isfield(S,'dZ')
            dZ = S.dZ;
        end
    end
end

if size(Out,2)*size(In,2) ~= size(Z,2)
    cond = size(Out,2)*size(In,2) == size(Z,2);
    cond = cond || size(Out,2)*size(In,2) == size(Z,2) - 1;
    assert(cond, 'size(Out,2) must equal size(Z,2)/size(In,2)');
end

Predicted = nan(size(Out));
Error = nan(size(Out));

logmsg('Computing Out_.Predicted, Out_.Error, PE, MSE, CC, SE, and Coherences.\n')
for k = 1:size(In,3) % Third dimension is segment
    for j = 1:size(Out,2) % Second dimension is component of Out

        u = k;
        if size(Z,3) == 1
            u = 1;
        end
        [Zi,~] = zinterp(fe,Z(:,:,u),size(In,1));
        if ~isempty(dZ)
            [dZi,~] = zinterp(fe,dZ(:,:,u),size(In,1));
        end

        zcols = (1:size(In,2)) + (j-1)*size(In,2);

        if ~isempty(dZ)
            Predicted(:,j,k) = zpredict(Zi(:,zcols),In(:,:,k),dZi);
        else
            Predicted(:,j,k) = zpredict(Zi(:,zcols),In(:,:,k));
        end

        Metrics.PE(1,j,k)  = pe_nonflag(Out(:,j,k), Predicted(:,j,k));
        Metrics.MSE(1,j,k) = mse_nonflag(Out(:,j,k), Predicted(:,j,k));
        Metrics.CC(1,j,k)  = cc_nonflag(Out(:,j,k), Predicted(:,j,k));

        [Metrics.PredictionCoherence(:,j,k), Metrics.fe] = ...
            coherence(Out(:,j,k), Predicted(:,j,k), 1, opts);

        Error(:,j,k) = Predicted(:,j,k)-Out(:,j,k);

        [Metrics.SN(:,j,k),~,Metrics.SNCLl(:,j,k),Metrics.SNCLu(:,j,k)] = ...
            signaltoerror(Out(:,j,k), Error(:,j,k), 1, opts);

        for i = 1:size(In,2)
            % If In and out both have two components, columns are
            %   <Outx,Inx>, <Outx,Iny>, <Outy,Iny>, <Outy,Iny>
            u = i + (j-1)*size(Out,2);
            Metrics.CrossCoherence(:,u,k) = ...
                coherence(Out(:,j,k), In(:,i,k), 1, opts);
        end

    end
end

logmsg('Computing DFT of Error and Predicted time series.\n')
for k = 1:size(In,3)
    DFTError(:,:,k) = dftbands(Error(:,:,k), opts);
    DFTPredicted(:,:,k) = dftbands(Predicted(:,:,k), opts);
end

logmsg('Computing residuals.\n');
for k = 1:size(In,3) % Third dimension is segment

    Zi = Z;
    dZi = dZ;

    if length(fe) ~= length(DFT.fe)
        logmsg('length(fe) ~= length(DFT.fe). Interpolating DFT onto fe\n');
        interpolate = 1;
    elseif ~all(fe == DFT.fe)
        logmsg('~all(fe == DFT.fe). Interpolating DFT onto fe\n');
        interpolate = 1;
    else
        interpolate = 0;
    end

    if interpolate
        Zi = zinterp(fe,Z,DFT.fe,{'linear', NaN});
        if ~isempty(dZi)
            dZi = zinterp(fe,dZ,DFT.fe,{'linear', NaN});
        end
    end

    for j = 1:size(Out,2) % Second dimension is component
        zcols = (1:size(In,2)) + (j-1)*size(In,2);
        for i = 1:length(DFT.In) % First dimension is frequency

            ftE = DFT.Out{i}(:,j,k);
            ftB = DFT.In{i}(:,:,k);
            Zf = Zi(i,zcols);
            if ~isempty(dZi)
                Zf = [Zf, dZi(i,j)];
                ftB = [ftB, ones(size(ftB,1),1)];
            end
            Metrics.Residuals{i,1}(:,j,k) =  ftE - ftB*transpose(Zf);
        end
    end

    if isfield(S,'Regression') && isfield(S.Regression,'Residuals')
        if size(S.Regression.Residuals,1) == size(Metrics.Residuals{i},1)
            logmsg('Checking residuals calculation.')
            for i = 1:length(S.Regression.Residuals)
                d = S.Regression.Residuals{i,:} - Metrics.Residuals{i,:};
                m(i,:) = max(abs(d(:)));
            end
            msg = 'Manually computed residuals do not match what what was returned by regression function.';
            assert(max(m) < 1e-10,msg);
        end
    end
end

if onsegments
    S.Segment.Metrics = Metrics;
    S.Segment.Out_.Error = Error;
    S.Segment.Out_.Predicted = Predicted;
    S.Segment.DFT.Out_.Error = DFTError;
    S.Segment.DFT.Out_.Predicted = DFTPredicted;
else
    S.Metrics = Metrics;
    S.Out_.Error = Error;
    S.Out_.Predicted = Predicted;
    S.DFT.Out_.Error = DFTError;
    S.DFT.Out_.Predicted = DFTPredicted;
end

if onsegments == 0 && isfield(S,'Segment')
    S = tflab_metrics(S,1);
end