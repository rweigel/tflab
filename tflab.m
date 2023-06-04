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
%  See also TRANSFERFNFD_OPTIONS, TRANSFERFNFD_TEST, TRANSFERFNFD_DEMO.

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
    return;
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

if isnan(opts.td.window.shift)
    % TODO: Allow window.width to be non-NaN?
    % No segmenting
    S = struct('In',B,'Out',E,'Options',opts);
    S = tflab_tdpreprocess(S);
    S = tflab_fdpreprocess(S);
    S = tflab_miso(S);
    S = tflab_metrics(S);
    return
end

% Compute TF for each segment of length opts.td.window.width
Tw = opts.td.window.width;
Ts = opts.td.window.shift;

assert(Tw > 1, 'opts.td.window.width must be greater than 1');
assert(Ts > 1, 'opts.td.window.shift must be greater than 1');

assert(Tw <= size(B,1), 'opts.td.window.width must be <= size(B,1)');
assert(mod(Tw,1) == 0, 'opts.td.window.width must be an integer');
assert(mod(Ts,1) == 0, 'opts.td.window.shift must be an integer');        

a = 1:Ts:size(B,1); % Start indices
b = a + Tw - 1;     % Stop indices

if b(end) > size(B,1)
    % If last window not full, omit it.
    a = a(1:end-1);
    b = b(1:end-1);
end
if b(end) ~= size(B,1)
    warning(['opts.td.window.width = %d,'...
             'opts.td.window.shift = %d, size(B,1) = %d.'...
             '\n\t Last %d point(s) will not be used.\n'],...
             Tw, Ts, size(B,1), size(B,1)-b(end));
end

% Compute TF for each segment
for s = 1:length(a)

    Iseg = a(s):b(s);
    if opts.tflab.loglevel > 0
        logmsg('Starting computations on segment %d of %d\n',s,length(a));
        logmsg('Segment time index range = [%d:%d]\n',Iseg(1),Iseg(end));
    end
    
    % Ss = Segment struct.
    Ss = struct();
    Ss.Options = opts;
    Ss.IndexRange = [a(s),b(s)]';

    Ss.In  = B(Iseg,:);
    Ss.Out = E(Iseg);

    Ss = tflab_tdpreprocess(Ss);
    Ss = tflab_fdpreprocess(Ss);

    if ~isempty(opts.fd.stack.average.function)
        % Only compute Z if performing stack averaging.

        logmsg('Computing Z for segment %d of %d\n',s,length(a));

        Ss = tflab_miso(Ss);

        if opts.tflab.loglevel > 0
            logmsg('Computing metrics on segment using segment Z.\n');
        end
        Ss = tflab_metrics(Ss);
        if opts.tflab.loglevel > 0
            logmsg('Computed metrics on segment using segment Z.\n');
        end

        if opts.tflab.loglevel > 0 && ~isempty(opts.fd.stack.average.function)
            logmsg('Computated Z and metrics for segment %d of %d; PE/CC/MSE %.2f/%.2f/%.3f\n',...
                   s,length(a),Ss.Metrics.PE,Ss.Metrics.CC,Ss.Metrics.MSE);
        end
    end
    
    if s == 1
        % First segment
        S = Ss;
    else
        % Combine segments
        S = combineStructs(S,Ss,3);
    end
end

if ~isempty(opts.fd.stack.average.function)
    if size(S.In,3) == 1
        S = rmfield(S,'IndexRange');
        % Only enough data for one segment, so no avg. needed.
        return;
    end
else
    % TF in a given frequency band is computed by doing regression on
    % DFTs in that frequency band for all segments.
    if isfield(opts.tflab, 'no_stack_regression')
        % Special case when E and B are cell arrays (non-continuous intervals).
        if opts.tflab.loglevel > 0
            logmsg(['opts.tflab.no_stack_regression set. '...
                    'Not doing stack regression (yet).\n']);
        end
        S.Segment = S;
        S = rmfield(S,'DFT');
        return;
    end    
end

% Move top-level structures under Segment
S.Segment = S;

% Remove Options from Segment b/c it applies to top-level
S.Segment = rmfield(S.Segment,'Options');

% Remove fields from top level of S b/c they were moved to S.Segments
fns = fieldnames(S);
for i = 1:length(fns)
    if ~strcmp('Segment',fns{i}) &&  ~strcmp('Options',fns{i})
        S = rmfield(S,fns{i});
    end
end

% Insert top-level structure elements
S.In = B;
S.Out = E;

if ~isempty(opts.fd.stack.average.function)

    % Compute TF by averaging segment TFs

    % Remove In and Out from Segment b/c they can be computed from IndexRange.
    S.Segment = rmfield(S.Segment,'In');
    S.Segment = rmfield(S.Segment,'Out');

    logmsg('Computing average of segment Zs.\n');
    S = stackAverage(S,opts);

    if opts.tflab.loglevel > 0
        logmsg('Computing metrics on unsegmented data using average Z.\n');
    end
    S = tflab_metrics(S);
    if opts.tflab.loglevel > 0
        logmsg('Computed metrics on unsegmented data using average Z.\n');
    end

else
    
    logmsg('Computing Z using stack regression.\n');        
    S = stackRegression(S, opts);

    % Remove In and Out from Segment b/c they can be computed from IndexRange.
    % (Then are needed for computing metrics in stackRegression, so must
    % be removed after stackRegression() call.)
    S.Segment = rmfield(S.Segment,'In');
    S.Segment = rmfield(S.Segment,'Out');

    logmsg('Computing metrics for full time series using computed Z.\n');
    S = tflab_metrics(S);
end

S.Options = opts;

end % tflab()

function S = stackAverage(S,opts)

    S.Z = mean(S.Segment.Z,3);
    S.fe = S.Segment.fe;

end

function S = stackRegression(S,opts)

    if opts.tflab.loglevel > 0
        logmsg('Starting stack regression.\n');
    end

    for i = 1:length(S.Segment.fe)

        tmp = squeeze(S.Segment.DFT.In(i,1,:));
        S.DFT.In{i,1} = cat(1,tmp{:});

        tmp = squeeze(S.Segment.DFT.Out(i,1,:));
        S.DFT.Out{i,1} = cat(1,tmp{:});

        tmp = squeeze(S.Segment.DFT.In(i,1,:));
        S.DFT.f{i,1} = cat(1,tmp{:});

        tmp = squeeze(S.Segment.DFT.Weights(i,1,:));
        S.DFT.Weights{i,1} = cat(1,tmp{:});
    end
    
    S.fe = S.Segment.fe;
    S = tflab_miso(S);
    if opts.tflab.loglevel > 0
        logmsg('Finished stack regression.\n');
    end
   
    % Temporarily insert Z computed using stack regression into Segment
    % because it will be used to compute metrics for each segment.
    S.Segment.Z = S.Z;

    if opts.tflab.loglevel > 0
        logmsg('Computing metrics on segments.\n');
    end
    
    S.Segment = tflab_metrics(S.Segment);
    
    % Remove inserted Z.
    S.Segment = rmfield(S.Segment,'Z');
end

function S = combineStructs(S1,S2,dim)
%combineStructs Combine tflab structures

    S = struct();    
    fns = fieldnames(S1);
    for i = 1:length(fns)
        if strcmp(fns{i},'Options')
            % Does not contain anything that can be concatenated and
            % if field exists, is same for S1 and S2.
            continue;
        end
        if isstruct(S1.(fns{i}))
            S.(fns{i}) = combineStructs(S1.(fns{i}),S2.(fns{i}),dim);
        else
            if (dim == 2 && strcmp(fns{i},'In')) || strcmp(fns{i},'fe')
                S.(fns{i}) = S1.(fns{i});
            else
                S.(fns{i}) = cat(dim,S1.(fns{i}),S2.(fns{i}));
            end
        end
    end
    % Use first structure's Options field. It will be same for S1 and S2.
    % If Options field does not exist, not at top-level.
    if isfield(S1,'Options')
        S.Options = S1.Options;
    end

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
    for c = 2:length(B)
        assert(size(B{c},2) == sB,...
            'Number of columns in B{%d} must be the same as B{1}.',c);
        assert(size(E{c},2) == sE,...
            'Number of columns in E{%d} must be the same as E{1}.',c);
        Nr(c) = size(B,1);
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
    S.Segment = struct();
    
    for c = 1:length(B) % Number of segments
        if opts.tflab.loglevel > 0
            logmsg('Starting computation for interval c = %d of %d\n',...
                    c,length(B));
        end
        % Compute Z for each segment in each interval
        if opts.tflab.loglevel > 0
            fprintf('%s\n',repmat('-',1,80));
            logmsg('Calling tflab(B{%d},E{%d},...)\n',c,c);
        end
        
        opts.tflab.no_stack_regression = 1;
        Sc = tflab(B{c},E{c},opts);
        
        Sc = rmfield(Sc,'In');
        Sc = rmfield(Sc,'Out');
        if ~isfield(Sc,'Segment')
            % If an interval had only one segment
            Sc.Segment = Sc;
        end
        
        if c == 1
            S = Sc;
        else
            % Combine segment fields across third dimension.
            S.Segment = combineStructs(S.Segment,Sc.Segment,3);
        end
    end
    
    if isempty(opts.fd.stack.average.function)
        
        % Remove fields from top level of S b/c they are calculations for
        % the last segment.
        fns = fieldnames(S);
        for i = 1:length(fns)
            if ~strcmp('Segment',fns{i})
                S = rmfield(S,fns{i});
            end
        end

        % Compute Z at each fe by regressing on all segment DFTs near fe.
        S = stackRegression(S, opts);

        % Calculate predicted/metrics/psd for each cell element        
        Sc = S;
        for c = 1:length(B)
            Sc.In = B{c};
            Sc.Out = E{c};
            
            if opts.tflab.loglevel > 0
                logmsg(...
                    ['Computing metrics for B{%d} and E{%d} using stack '...
                     'regression transfer function.\n'],c,c);
            end
            Sc = tflab_metrics(Sc);
            if opts.tflab.loglevel > 0
                logmsg(...
                    ['Finished computing metrics for B{%d} and E{%d} using '...
                     'stack regression transfer function.\n'],c,c);
            end
            for j = 1:size(Sc.Metrics.PE, 2)
                logmsg( ...
                        'Output col %d, PE/CC/MSE = %.2f/%.2f/%.3f\n',...
                         j,...
                         Sc.Metrics.PE(1,j),...
                         Sc.Metrics.CC(1,j),...
                         Sc.Metrics.MSE(1,j));
            end
            
            S.Metrics{c} = Sc.Metrics;
        end
        
        S.In = B;
        S.Out = E;
    else       
        % Compute stack average Z and its predicted/metrics/psd for
        % full In/Out
        S.In = B;
        S.Out = E;
        S = rmfield(S,'DFT');
        S = rmfield(S,'Regression');
        S = rmfield(S,'Z');
        logmsg('Computing stack average Z and its metrics for full In/Out.\n');
        S = stackAverage(S,opts);
    end
end