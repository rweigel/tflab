function S = tflab(B,E,t,opts)
%TRANSFERFNFD Frequency domain MIMO transfer function estimate
%
%  S = tflab(B,E) returns a structure with an estimate of the
%  transfer function Z in the expression
%
%    Ex(f) = Zxx(f)Bx(f) + Zxy(f)By(f) + ...
%  
%  given time series for Ex(t), Bx(t), By(t), ... and using the convention
%  that for an arbitrary variable U, U(f) is the fourier transform of U(t).
%
%  Estimates are made for the complex-valued transfer function Z at a set
%  of evaluation frequencies, fe, using regression with a set of
%  frequencies in a band around each fe on the model equation above.
%  
%  The set of evaluation frequencies and windows are determined using the
%  function evalfreq(). By default, the evaluation frequencies are
%  lograrithmically spaced with approximately 7 frequencies per decade.
%
%  If the number of columns in B is Nb and E has Ne columns, Z will have
%  2*Nb*Ne columns with columns 1:Nb corresponding to the transfer function
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
%  S = tflab(B,E,t) associates a time value with each row of B and
%  E. Use t = [] to use the default of t = [1:size(B,1)]'.
%
%  S = tflab(B,E,t,options) uses options returned by the function
%  tflab_options.
%
%  See also TRANSFERFNFD_OPTIONS, TRANSFERFNFD_TEST, TRANSFERFNFD_DEMO.

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

if nargin == 2
    t = [];
    opts = [];
end

if nargin == 3
    % tflab(B,E,t) or
    % tflab(B,E,opts)
    if isstruct(t)
        % tflab(B,E,opts)
        opts = t;
        t = [];
    else
        % tflab(B,E,t)        
        opts = [];
        assert(size(t,1) == size(B,1),...
                'Required: size(t,1) == size(B,1) == size(E,1)');    
    end
end

if nargin < 3 || isempty(opts)
    % tflab(B,E) or
    % tflab(B,E,t)
    % Use default options.
    opts = tflab_options(1);
    if opts.tflab.loglevel
        logmsg(['No options given. '...
                 'Using options returned by tflab_options(1)\n']);
    end
end

% For each matrix in B (and corresponding matrix in E), do TF calculations.
if iscell(B)
    
    % Each cell contians an interval and an arbitrary gap in time is
    % assumed between the end of one cell and the start of the next
    % cell.
    assert(all(size(B) == size(E)),'Required: size(B) == size(E)');
    assert(isvector(B),'Required: B must be vector cell array');
    assert(isvector(E),'Required: E must be vector cell array');
    if ~isempty(t)
        assert(all(size(t) == size(B)),'Required: size(t) == size(B)');
        assert(isvector(t),'Required: t must be vector cell array');
    else
        for c = 1:length(B)
            t{c} = (opts.td.start:size(B{c},1))';
        end
    end
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
            logmsg(...
                'Starting computation for disconnected segment c = %d of %d\n',...
                c,length(B));
        end
        % Compute Z for each segment in each interval
        if opts.tflab.loglevel > 0
            fprintf('%s\n',repmat('-',1,80));
            logmsg('Calling tflab(B{%d},E{%d},...)\n',c,c);
        end
        
        opts.tflab.no_stack_regression = 1;
        Sc = tflab(B{c},E{c},t{c},opts);
        
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
        %S.Segment = stackRegression(S.Segment, opts);
        S = stackRegression(S, opts);
                
        % Calculate predicted/metrics/psd for each cell element        
        Sc = S;
        for c = 1:length(B)
            Sc.In = B{c};
            Sc.Out = E{c};
            Sc.Time = t{c};

            if opts.tflab.loglevel > 0
                logmsg(...
                    ['Computing metrics for B{%d} and E{%d} using stack '...
                     'regression transfer function.\n'],c,c);
            end
            Sc = tflab_metrics(Sc,opts);
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
        S.Time = t;
    else       
        % Compute stack average Z and its predicted/metrics/psd for
        % full In/Out
        S.In = B;
        S.Out = E;
        S.Time = t;
        S = rmfield(S,'DFT');
        S = rmfield(S,'Regression');
        S = rmfield(S,'Z');
        logmsg('Computing stack average Z and its metrics for full In/Out.\n');
        S = stackAverage(S,opts);
    end
    return
end

% Number of time values must be the same.
assert(size(B,1) == size(E,1),...
        'Required: size(B,1) == size(E,1)');

assert(size(B,1) >= size(B,2),...
    ['Not enough time samples: size(B,1) must be greater than '...
     'or equal to size(B,2)']);

if nargin < 3 || isempty(t)
    t = (opts.td.start:size(B,1))';
end

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
        Sc = tflab(B,E(:,j),t,opts);
        if j == 1
            S = Sc;
        else
            S = combineStructs(S,Sc,2);
        end
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code start. Given Nt x Nin B and Ntx1 E, compute TF, where Nt is the
% number of timesteps and Nin is the number of inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.tflab.loglevel > 0
    logmsg( ['Computing transfer function for input/output '...
                    'sizes [%d,%d]/[%d,1]\n'],...
                     size(B),size(E,1));
    if opts.tflab.loglevel > 1
        logmsg( 'Options:\n');
        printstruct(opts);
    end
end

if isnan(opts.td.window.width)
    % Compute one TF (no segmenting)
    if opts.td.window.loglevel
        logmsg('opts.td.window.with is NaN. Using size(B,1).\n');
    end
    % Set default window width and shift equal to the number of time points.
    opts.td.window.width = size(B,1);
    opts.td.window.shift = size(B,1);
    S = tflab(B,E,t,opts);
    return
end

% Compute TF for each segment of length opts.td.window.width
assert(opts.td.window.width > 1,...
        'opts.td.window.width must be greater than 1');
assert(opts.td.window.shift > 1,...
        'opts.td.window.shift must be greater than 1');
assert(opts.td.window.width <= size(B,1),...
        'opts.td.window.width must be less than or equal to size(B,1)');
Tw = opts.td.window.width;
Ts = opts.td.window.shift;
assert(mod(Tw,1) == 0,...
        'opts.td.window.width must be an integer');
assert(mod(Ts,1) == 0,...
        'opts.td.window.shift must be an integer');        
a = 1:Ts:size(B,1); % Start indices
b = a + Tw - 1;     % Stop indices
if b(end) > size(B,1)
    % If last window not full, omit it.
    a = a(1:end-1);
    b = b(1:end-1);
end
if b(end) ~= size(B,1)
    warning(...
        ['opts.td.window.width = %d,'...
         'opts.td.window.shift = %d, size(B,1) = %d.'...
         '\n\t Last %d point(s) will not be used.\n'],...
         Tw,Ts,size(B,1),size(B,1)-b(end));
end
optsx = opts;
optsx.td.window.width = NaN;
optsx.td.window.shift = NaN;
% Compute TF for each segment
for s = 1:length(a)
    Iseg = a(s):b(s);
    if opts.tflab.loglevel > 0
        logmsg(...
                'Starting computation for segment %d of %d\n',...
                s,length(a));
    end
    % Ss = Segment struct.
    %Ss = tflab(B(Iseg,:),E(Iseg,:),t(Iseg),optsx);
    Ss = main(B(Iseg,:),E(Iseg,:),t(Iseg),optsx);
    Ss.IndexRange = [a(s),b(s)]';
    if opts.tflab.loglevel > 0 ...
                    && ~isempty(opts.fd.stack.average.function)
        % Summarize results for each column of E
        for j = 1:size(E,2)
            logmsg(...
                    ['Finished segment %d of %d; PE/CC/MSE '...
                     'of Out(%d:%d,%d) = %.2f/%.2f/%.3f\n'],...
                     s,...
                     length(a),...
                     Iseg(1),...
                     Iseg(end),...
                     j,...
                     Ss.Metrics.PE(j),...
                     Ss.Metrics.CC(j),...
                     Ss.Metrics.MSE(j));
        end
    end
    if s == 1
        S = Ss;
    else
        S = combineStructs(S,Ss,3);
    end
end

if ~isempty(opts.fd.stack.average.function)
    % TF is average of segment TFs
    if size(S.In,3) > 1 % More than one segment
        % Move top-level structures under Segment
        S.Segment = S;
        S.Segment = rmfield(S.Segment,'Options');
        S = rmfield(S,'Metrics');
        if isfield(S,'Window')
            S = rmfield(S,'Window');
        end
        S = rmfield(S,'IndexRange');
        S = rmfield(S,'DFT');
        S = rmfield(S,'Regression');
        S = rmfield(S,'Z');
        S = rmfield(S,'Phi');
        S = rmfield(S,'H');
        S = rmfield(S,'tH');
        % Insert top-level structure elements
        S.In = B;
        S.Out = E;
        S.Time = t;
        logmsg('Computing stack average Z.\n');
        S = stackAverage(S,opts);
    end
else
    % TF in a given frequency band is computed by doing regression on
    % DFTs in that frequency band for all segments.
    if isfield(opts.tflab, 'no_stack_regression')
        if opts.tflab.loglevel > 0
            logmsg(...
                ['opts.tflab.no_stack_regression set. '...
                 'Not doing stack regression (yet).\n']);
        end
        S.Segment = S;
        S = rmfield(S,'DFT');
        return
    else

        % Compute Z using stack regression. Computes segment metrics
        % using computed Z.
        logmsg('Computing Z using stack regression.\n');
        S.Segment = S;
        S.Segment = rmfield(S.Segment,'Options');

        % Remove fields from top level of S b/c they are calculations for
        % all segments.
        fns = fieldnames(S);
        for i = 1:length(fns)
            if ~strcmp('Segment',fns{i}) &&  ~strcmp('Options',fns{i})
                S = rmfield(S,fns{i});
            end
        end
        S = stackRegression(S, opts);

        % Insert top-level structure elements
        S.In = B;
        S.Out = E;
        S.Time = t;

        logmsg('Computing stack regression metrics.\n');
        S = tflab_metrics(S,opts);
    end
end

S.Options = opts;

function S = main(B, E, t, opts)

    S = struct();
        S.In = B;
        S.Out = E;
        S.Time = t;
        S.Options = opts;
        if ~isempty(opts.fd.stack.average.function)
            S.Regression = struct();
        end

    if ~isempty(opts.td.detrend.function)
        B = opts.td.detrend.function(B,opts.td.detrend.functionargs{:});
        E = opts.td.detrend.function(E,opts.td.detrend.functionargs{:});
    end
    
    % TODO?: Allow TD window and prewhiten to not be same for input and output
    % and then compute corrected Z.
    if ~isempty(opts.td.window.function)
        if opts.tflab.loglevel > 0
            logmsg('Windowing input and output using %s\n',opts.td.window.functionstr);
        end
        [B,~] = opts.td.window.function(B,opts.td.window.functionargs{:});
        [E,W] = opts.td.window.function(E,opts.td.window.functionargs{:});
        S.Window.Weights = W;
        S.Window.In = B;
        S.Window.Out = E;    
    else
        if opts.tflab.loglevel > 0
            logmsg( ...
                'No time domain window applied b/c no function given.\n');
        end
    end

    if ~isempty(opts.td.prewhiten.function)
        if opts.tflab.loglevel > 0
            logmsg('Prewhitening input and output using %s\n',opts.td.prewhiten.functionstr);
        end
        S.Prewhiten = struct();
        S.Prewhiten.Comment = 'S.Prewhiten.In = prewhitening applied to S.Window.In (or S.In if now time domain window given).';

        [B,a,b] = opts.td.prewhiten.function(B,opts.td.prewhiten.functionargs{:});
        S.Prewhiten.InFilter = [a,b];
        S.Prewhiten.In = B;

        [E,a,b] = opts.td.prewhiten.function(E,opts.td.prewhiten.functionargs{:});
        S.Prewhiten.OutFilter = [a,b];
        S.Prewhiten.Out = E;
    else
        if opts.td.prewhiten.loglevel
            logmsg( ...
                    ['No time domain prewhitening performed '...
                     'b/c no function given.\n']);
        end
    end

    if ~isnan(opts.td.zeropad)
        if opts.tflab.loglevel > 0
            logmsg('Zero padding input and output with %d zeros\n',opts.td.zeropad);
        end
        E = [E;zeros(opts.td.zeropad,1)];
        B = [B;zeros(opts.td.zeropad,size(B,2))];
        S.Zeropad.In = E;
        S.Zeropad.Out = B;
        S.Zeropad.N = opts.td.zeropad;
    end

    [~,f] = fftfreq(size(B,1)); % Unique DFT frequencies

    if opts.tflab.loglevel > 0
        logmsg(['Calling %s() using additional arguments given in \n'...
                'opts.fd.evalfreq.functionargs\n'],...
                func2str(opts.fd.evalfreq.function));
    end
    [fe,Ic,Ne] = opts.fd.evalfreq.function(...
                    size(B,1),opts.fd.evalfreq.functionargs{:});
    S.fe = fe';


    if opts.tflab.loglevel > 0
        logmsg( 'Computing raw DFTs of input and output.\n');
    end

    ftB = fft(B);
    ftE = fft(E);

    % Compute # of unique frequency values.
    N = size(B,1);
    if mod(N,2) == 0
        Np = N/2 + 1; % f = -0.5 value is kept.
    else
        Np = (N-1)/2 + 1;    
    end

    ftB = ftB(1:Np,:);
    ftE = ftE(1:Np,:);

    if opts.tflab.loglevel > 0
        if isempty(opts.fd.stack.average.function) || ~strcmp(opts.fd.program.name,'tflab')
            logmsg(...
                'Starting freq band calcs for %d frequencies.\n',...
                length(Ic)-1);
            logmsg(...
                ['opts.fd.stack.average.function = '''' =>\n' ...
                'No regression performed for each freq. band of segment\n']);
        else
            logmsg(['Starting freq band and regression '...
                            'calcs for %d frequencies.\n'],length(Ic));
            logmsg(...
                ['Using %s() with additional arguments given in\n'...
                 'opts.fd.regression.functionargs\n'],...
                func2str(opts.fd.regression.function));
        end
    end

    winfn = opts.fd.window.function;
    if opts.fd.window.loglevel && strcmp(opts.fd.program.name,'tflab')
        logmsg( 'Using FD window function %s\n',func2str(winfn));
    end

    for j = 1:length(Ic)

        if opts.fd.regression.loglevel ...
                && ~isempty(opts.fd.stack.average.function) ...
                && strcmp(opts.fd.program.name,'tflab')
            logmsg(...
                    ['Starting freq band and regression '...
                     'calcs on frequency %d of %d\n'],...
                     j, length(fe));
        end

        W = winfn(2*Ne(j)+1);
        W = W/sum(W);
        r = Ic(j)-Ne(j):Ic(j)+Ne(j); % Index range

        W  = sqrt(W);
        Wr = repmat(W,1,size(ftB,2));

        S.DFT.Out{j,1} = ftE(r,1);
        S.DFT.In{j,1} = ftB(r,:);
        S.DFT.f{j,1} = f(r);
        S.DFT.Weights{j,1} = W;

        if opts.fd.window.loglevel
            logmsg(...
                    ['Band with center of fe = %.8f has %d '...
                     'points; fl = %.8f fh = %.8f\n'],...
                     fe(j),...
                     length(r),...
                     f(Ic(j)-Ne(j)),...
                     f(Ic(j)+Ne(j)));
        end

        % If not computing Z based on stack averages, don't need to do
        % regression as it is done later.
        if ~isempty(opts.fd.stack.average.function) ...
                && strcmp(opts.fd.program.name,'tflab')
            %args = opts.fd.regression.functionargs;    
            regressargs = opts.fd.regression.functionargs;
            regressfunc = opts.fd.regression.function;
        
            lastwarn('');

            if length(r) < size(ftB,2)
                if fe(j) == 0
                    Z(j,:) = zeros(1,size(ftB,2));
                    S.Regression.Weights{j,1} = nan*W;
                    S.Regression.Residuals{j,1} = nan*W;
                    logmsg('! System if underdetermined for fe = %f. Setting Z equal to zero(s) for this frequency.',fe(j));
                    %warning(sprintf('System if underdetermined for fe = %f. Setting Z equal to zero(s) for this frequency.',fe(j)));
                else
                    Z(j,:) = nan(1,size(ftB,2));
                    S.Regression.Weights{j,1} = nan*W;
                    S.Regression.Residuals{j,1} = nan*W;
                    logmsg('! System if underdetermined for fe = %f. Setting Z equal to NaN(s) for this frequency.',fe(j));
                    %warning(sprintf('System if underdetermined for fe = %f. Setting Z equal to NaN(s) for this frequency.',fe(j)));
                end
                continue;
            end

            if length(r) == 1 && size(ftB,2) == 1
                if ftB(r,1) == 0 && ftE(r,1) ~= 0 && fe(j) == 0.5
                    % Special case for when there is an evaluation freq.
                    % at 0.5. Will occur if frequency window is of length
                    % 1, as it is for some of the tests and demos.
                    Z(j,:) = zeros(1,size(ftB,2)) + 1j*zeros(1,size(ftB,2));
                    S.Regression.Weights{j,1} = nan*W;
                    S.Regression.Residuals{j,1} = nan*W;
                    logmsg('! System if underdetermined for fe = %f. Setting Z equal to zero(s) for this frequency.',fe(j));
                    continue;
                end
            end

            %[Z(j,:),Weights,Residuals] = ...
            %                opts.fd.regression.function(...
            %                      Wr.*ftB(r,:),W.*ftE(r,1),args{:});
            [Z(j,:),Residuals,Weights] = ...
                       regressfunc(W.*ftE(r,1),Wr.*ftB(r,:),regressargs{:});

            if ~isempty(lastwarn)
                logmsg('Above is for eval. freq. #%d; fe = %f; Te = %f\n', ...
                    j,fe(j),1/fe(j));
                logmsg(sprintf('ftE = \n'));
                logmsg(sprintf('   %.16f\n',ftE(r,1)));
                logmsg(sprintf('ftB = \n'));
                logmsg(sprintf('   %.16f\n',ftB(r,:)));
            end
            
            S.Regression.Weights{j,1} = Weights;
            S.Regression.Residuals{j,1} = Residuals;
        end

    end

    if opts.tflab.loglevel > 0
        if opts.fd.regression.loglevel ...
                && ~isempty(opts.fd.stack.average.function) ...
                && strcmp(opts.fd.program.name,'tflab')
            logmsg(...
                ['Finished freq band and regression '...
                 'calculations for %d eval. freqs.\n'],...
                 length(Ic)-1);

        end
    end

    if isempty(opts.fd.stack.average.function)
        % If a stack average function was not given, then a Z for each
        % segment was not computed.

    else
        % If a stack average function is given, then a Z for each segment
        % is computed and Z is the stack average of the segment Zs.

        if all(isnan(Z(:)))
            error('All Z values are NaN');
        end

        S.Z = Z;

        if opts.tflab.loglevel > 0
            logmsg('Computing Phi\n');
        end
        S.Phi = atan2(imag(S.Z),real(S.Z));
        if opts.tflab.loglevel > 0
            logmsg('Computed Phi\n');
        end

        if opts.tflab.loglevel > 0
            logmsg('Interpolating Z\n');
        end
        [Zi,~,Zir,fir] = zinterp(S.fe,S.Z,size(S.In,1));
        %S.Zi = Zir;
        %S.fi = fir;
        if opts.tflab.loglevel > 0
            logmsg('Interpolated Z\n');
        end

        if opts.tflab.loglevel > 0    
            logmsg('Computing H\n');
        end
        [S.H,S.tH] = z2h(Zi);
        if opts.tflab.loglevel > 0    
            logmsg('Computed H\n');
        end

        if opts.tflab.loglevel > 0
            logmsg('Computing metrics\n');
        end
        S = tflab_metrics(S,opts);
        if opts.tflab.loglevel > 0
            logmsg('Computed metrics\n');
            logmsg('PE/CC/MSE = %.2f/%.2f/%.3f\n',...
                     S.Metrics.PE,...
                     S.Metrics.CC,...
                     S.Metrics.MSE);
        end

    end
end % main()
end % tflab()

function S = stackAverage(S,opts)

    S.Z = mean(S.Segment.Z,3);
    S.Phi = atan2(imag(S.Z),real(S.Z));
    
    if iscell(S.In)
        S.Metrics = {};
        for i = 1:length(S.In)
            tmp = struct();
            tmp.In = S.In{i};
            tmp.Out = S.Out{i};
            tmp.Z = S.Z;
            tmp.fe = S.fe;
            tmp = tflab_metrics(tmp,opts);
            S.Metrics{i} = tmp.Metrics;
        end
    else
        if opts.tflab.loglevel > 0
            logmsg('Computing metrics using stack averaged Z.\n');
        end        
        S = tflab_metrics(S,opts);
        if opts.tflab.loglevel > 0
            logmsg('Computed metrics using stack averaged Z.\n');
        end
    end
end

function S = stackRegression(S,opts)

    % At each evaluation frequency index i and output column c and for each
    % segment s, S.DFT.Out(i,c,s) is a single-column matrix with rows of
    % DFTs of segment s in the frequency band associated with i.
    % S.DFT.In(i,1,s) is a matrix with same number of columns of S.In. Each
    % column of S.DFT.In(i,1,s) contains the DFTs for the respective column
    % in S.In for freq. band i and segment s.

    if opts.tflab.loglevel > 0
        logmsg('Starting stack regression.\n');
    end

    regressargs = opts.fd.regression.functionargs;
    regressfunc = opts.fd.regression.function;
    
    for i = 1:size(S.Segment.DFT.In, 1) % Eval frequencies
        for c = 1:size(S.Segment.DFT.Out, 2) % Columns of E

            if opts.tflab.loglevel > 1
                logmsg('Performing stack regression for eval freq. %d and on column %d of input.\n', i, c);
            end
            
            tmp = squeeze(S.Segment.DFT.In(i,1,:));
            ftB = cat(1,tmp{:});

            tmp = squeeze(S.Segment.DFT.Out(i,c,:));
            ftE = cat(1,tmp{:});

            tmp = squeeze(S.Segment.DFT.Weights(i,c,:));
            W   = cat(1,tmp{:});
            Wr  = repmat(W,1,size(ftB,2));
            if length(W) < size(ftB,2)
                if S.Segment.fe(i) == 0
                    z = zeros(1,size(ftB,2));
                    warning(sprintf('System if underdetermined for fe = %f. Setting Z equal to zero(s) for this frequency.',S.Segment.fe(i)));
                else
                    z = nan(1,size(ftB,2));
                    warning(sprintf('System if underdetermined for fe = %f. Setting Z equal to NaN(s) for this frequency.',S.Segment.fe(i)));
                end
                S.Regression.Weights{i,c} = nan*W;
                S.Regression.Residuals{i,c} = nan*W;
            else          
                % https://www.mathworks.com/matlabcentral/answers/364719-detect-warning-and-take-action#answer_289064    
                warning('');
                %[z,weights,residuals] = opts.fd.regression.function(Wr.*ftB,W.*ftE,args{:});
                [z,residuals,weights] = regressfunc(W.*ftE,Wr.*ftB,regressargs{:});                warnMsg = lastwarn;
                if ~isempty(warnMsg)
                    logmsg('Warning above occured on output column %d, eval freq. = %g\n', c, S.Segment.fe(i));
                    %ftE
                    %ftB
                    %keyboard
                end
                S.DFT.f{i,c} = S.Segment.DFT.f{i,c};
                S.DFT.Out{i,c} = ftE;
                S.DFT.In{i,c} = ftB;
                S.Regression.Weights{i,c} = weights;
                S.Regression.Residuals{i,c} = residuals;
            end
            if c == 1
                Zc = z.';
            else
                Zc = [Zc,z.'];
            end
            if opts.tflab.loglevel > 1
                logmsg('Performed stack regression for eval freq. %d and on column %d of input.\n', i, c);
            end            
        end
        Z(i,:) = Zc;
    end

    if opts.tflab.loglevel > 0
        logmsg('Finished stack regression.\n');
    end
    
    S.fe = S.Segment.fe;
    S.Z = Z;    

    if opts.tflab.loglevel > 1
        logmsg('Computing Phi\n');
    end
    S.Phi = atan2(imag(Z),real(Z));
    if opts.tflab.loglevel > 1
        logmsg('Finished computing Phi\n');
    end

    if opts.tflab.loglevel > 1
        logmsg('Interpolating Z\n');    
    end
    [Zi,~,Zir,fir] = zinterp(S.fe,S.Z,size(S.Segment.In,1));
    %S.Zi = Zir;
    %S.fi = fir;
    if opts.tflab.loglevel > 1
        logmsg('Finished interpolating Z\n');
    end

    if opts.tflab.loglevel > 1
        logmsg('Computing H\n');
    end
    [S.H,S.tH] = z2h(Zi);
    if opts.tflab.loglevel > 1
        logmsg('Finished computing H\n');
    end

    if opts.tflab.loglevel > 0
        logmsg(...
                ['Computing metrics on each segment using stack regression '...
                 'transfer function.\n']);
    end
    
    S.Segment.Z = Z;
    S.Segment = tflab_metrics(S.Segment,opts);
    S.Segment = rmfield(S.Segment,'Z');

    if opts.tflab.loglevel > 0
        for c = 1:size(S.Segment.Metrics.PE, 2)
            for s = 1:size(S.Segment.Metrics.PE, 3)
                logmsg(...
                        'Segment %d: PE/CC/MSE = %.2f/%.2f/%.3f\n',...
                         s,...
                         S.Segment.Metrics.PE(1,c,s),...
                         S.Segment.Metrics.CC(1,c,s),...
                         S.Segment.Metrics.MSE(1,c,s));
            end
        end
        logmsg(...
                ['Finished computing metrics on each segment using stack regression '...
                 'transfer function.\n']);
    end

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