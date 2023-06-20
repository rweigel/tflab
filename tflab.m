function S = tflab(B,E,opts)
%TFLAB Wrapper function for tflab_miso.
%
%  Frequency domain MIMO transfer function estimate
%
%  tf = TFLAB(B,E)
%
%  If B has Nb and E has Ne columns, Z will have 2*Nb*Ne columns, with
%  columns 1:Nb corresponding to the transfer function
%  that gives E(:,1) given B, columns Nb+1:2*Nb corresponding to the
%  transfer function that gives E(:,2) given B, etc.
%
%  In more familar notation, if the input and outputs are
%
%   Bx(t) = B(:,1)  By(t) = B(:,2)
%   Ex(t) = E(:,1)  Ey(t) = E(:,2)
%
%  the model equations are
%
%   Ex(f) = Zxx(f)Bx(f) + Zxy(f)By(f)
%   Ey(f) = Zyx(f)Bx(f) + Zyy(f)By(f)
%
%  and S = tflab(B,E) returns a structure S with a field Z such that
%
%   Zxx = Z(:,1)   Zxy = Z(:,2)
%   Zyx = Z(:,3)   Zyy = Z(:,4)
% 
%  with the rows of Z being estimates of the transfer function at the
%  evaluation frequencies given by field fe in S.
%
%  S = tflab(B,E,options) uses options returned by the function
%  tflab_options.
%
%  See also TFLAB, TFLAB_TEST, TFLAB_DEMO.

addpath(fullfile(fileparts(mfilename('fullpath'))));
tflab_setpaths();

if nargin == 2
    % tflab(B,E) or
    % Use default options.
    opts = tflab_options(1);
    if opts.tflab.loglevel
        logmsg(['No options given. '...
                 'Using options returned by tflab_options(1)\n']);
    end
end

% For each matrix in B (and corresponding matrix in E), do TF calculations.
if iscell(B)
    S = tflab_intervals(B, E, opts);
    return
end

% Number of time values must be the same.
assert(size(B,1) == size(E,1),'Required: size(B,1) == size(E,1)');

assert(size(B,1) >= size(B,2),...
    ['Not enough time samples: size(B,1) must be greater than '...
     'or equal to size(B,2)']);

% If E has more than one column, do TF calculation for each column.
if size(E,2) > 1
    % Ex = ZxxBx + ZxyBy + ...
    % Ey = ZyxBx + ZyyBy + ...
    % ...
    for j = 1:size(E,2)
        if opts.tflab.loglevel > 0
            fprintf('%s\n',repmat('-',1,80));
            logmsg('Calling tflab(B,E(:,%d),...)\n',j);
        end
        Sc = tflab(B,E(:,j),opts);
        if j == 1
            S = Sc;
        else
            S = combineStructs(S,Sc,2);
        end
    end
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code start. Given size(B) = [Nt, Nin] and size(E) = [Nt,1], compute
% TF, where Nt is the number of timesteps and Nin is the number of inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(isnan(E))
    error('E has NaNs');
end
if any(isnan(B))
    error('B has NaNs');
end

if opts.tflab.loglevel > 0
    logmsg(['Computing transfer function for input/output '...
            'sizes [%d, %d]/[%d, 1]\n'],size(B),size(E,1));
    if opts.tflab.loglevel > 1
        logmsg('Options:\n');
        printstruct(opts);
    end
end

if isnan(opts.td.window.width) || size(E,1) == opts.td.window.width
    % No segmenting
    S = struct('In',B,'Out',E,'Options',opts);
    S = tflab_tdpreprocess(S);
    S = tflab_fdpreprocess(S);
    [S.Z,S.fe,S.Regression] = tflab_miso(S.DFT,opts);
    S = tflab_metrics(S);
    return
end

Sc = tflab_preprocess(B,E,opts);

S = struct('In',B,'Out',E,'Options',opts);

if ~isempty(opts.fd.stack.average.function)

    logmsg('Computing average of segment Zs.\n');
    [S.Z,S.fe,S.Segment] = stackAverage(Sc,opts);

    if opts.tflab.loglevel > 0
        logmsg('Computing metrics on unsegmented data using average Z.\n');
    end
    S = tflab_metrics(S);
    if opts.tflab.loglevel > 0
        logmsg('Computed metrics on unsegmented data using average Z.\n');
    end
else

    logmsg('Computing Z using stack regression.\n');        
    [S.Z,S.fe,S.Segment,S.Regression] = stackRegression(Sc, opts);

    logmsg('Computing metrics for full time series using computed Z.\n');
    S = tflab_metrics(S);
end

end % tflab()

function [Z,fe,Segment] = stackAverage(Sc,opts)

    for s = 1:length(Sc)

        if opts.tflab.loglevel > 0
            if isfield(Sc{s},'IndexRange')
                a = Sc{s}.IndexRange(1);
                b = Sc{s}.IndexRange(2);
            else
                a = 1;
                b = size(Sc{s}.In,1);
            end
            logmsg('Starting computations on segment %d of %d\n',s,length(Sc));
            logmsg('Segment time index range = [%d:%d]\n',a,b);
        end

        logmsg('Computing Z for segment %d of %d\n',s,length(Sc));
        [Sc{s}.Z,Sc{s}.fe,Sc{s}.Regression] = tflab_miso(Sc{s}.DFT,opts);

        if opts.tflab.loglevel > 0
            logmsg('Computing metrics on segment using segment Z.\n');
        end
        Sc{s} = tflab_metrics(Sc{s});
        if opts.tflab.loglevel > 0
            logmsg('Computed metrics on segment using segment Z.\n');
        end

        if opts.tflab.loglevel > 0 && ~isempty(opts.fd.stack.average.function)
            logmsg('Computated Z and metrics for segment %d of %d; PE/CC/MSE %.2f/%.2f/%.3f\n',...
                   s,length(Sc),Sc{s}.Metrics.PE,Sc{s}.Metrics.CC,Sc{s}.Metrics.MSE);
        end
    end    

    Segment = combineStructs(Sc,3);

    Z = mean(Segment.Z,3);
    fe = Segment.fe;

end

function [Z,fe,Segment,Regression] = stackRegression(Sc,opts)

    Segment = combineStructs(Sc,3);

    if opts.tflab.loglevel > 0
        logmsg('Starting stack regression.\n');
    end
    
    DFT = dftcombine(Segment.DFT);
    [S.Z,S.fe,S.Regression] = tflab_miso(DFT,opts);
    
    if opts.tflab.loglevel > 0
        logmsg('Finished stack regression.\n');
    end
   
    % Temporarily insert Z computed using stack regression into Segment
    % because it will be used to compute metrics for each segment.
    Segment.Z = S.Z;
    Segment.fe = S.fe;

    if opts.tflab.loglevel > 0
        logmsg('Computing metrics on segments.\n');
    end

    Segment = tflab_metrics(Segment);

    % Remove inserted Z.
    Segment = rmfield(Segment,'Z');
    Segment = rmfield(Segment,'fe');

    Z = S.Z;
    fe = S.fe;
    Regression = S.Regression;
    
end

function S = tflab_intervals(B, E, opts)    

    % Each cell contians an interval and an arbitrary gap in time is
    % assumed between the end of one cell and the start of the next
    % cell.

    assert(all(size(B) == size(E)),'Required: size(B) == size(E)');
    assert(isvector(B),'Required: B must be vector cell array');
    assert(isvector(E),'Required: E must be vector cell array');
    sE = size(E{1},2);
    sB = size(B{1},2);

    % Check that number of columns is same for all intervals.
    for i = 2:length(B)
        assert(size(B{i},2) == sB,...
            'Number of columns in B{%d} must be the same as B{1}.',i);
        assert(size(E{i},2) == sE,...
            'Number of columns in E{%d} must be the same as E{1}.',i);
        Nr(i) = size(B,1);
    end
    
    if isnan(opts.td.window.width) && length(unique(Nr)) ~= 1
        warning(['opts.td.window.width not given and intervals do not', ...
                'all have the same length.\n', ...
                'Using shortest interval length (%d) for width and shift'],...
                min(Nr));
        opts.td.window.width = min(Nr);
        opts.td.window.shift = min(Nr);
    end
    
    S = struct();
    
    c = 1;
    for i = 1:length(B) % Number of intervals
        Si = tflab_preprocess(B{i},E{i},opts);
        if isstruct(Si)
            Si = {Si};
        end
        for s = 1:length(Si)
            Si{s}.SegmentNumber = i;
            Sc{c} = Si{s};
            c = c+1;
        end
    end

    S.In = B;
    S.Out = E;
    S.Options = opts;

    if ~isempty(opts.fd.stack.average.function)
        logmsg('Computing average of segment Zs.\n');
        [S.Z,S.fe,S.Segment] = stackAverage(Sc,opts);
    end
    
    if isempty(opts.fd.stack.average.function)
        logmsg('Computing Z using stack regression.\n');        
        [S.Z,S.fe,S.Segment,S.Regression] = stackRegression(Sc, opts);
    end
end