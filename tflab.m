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
    S = intervals_(B, E, opts);
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
    return
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
    S = tflab_tdpreprocess(S,opts.tflab.loglevel);
    S = tflab_fdpreprocess(S);
    [S.Z,S.fe,S.dZ,S.Regression] = tflab_miso(S.DFT,opts);
    S = tflab_metrics(S);
    return
end

if ~isempty(opts.fd.stack.average.function)

    S = tflab_preprocess(B,E,opts,'both',1,0);
    S = stackaverage_(S,opts);
    
    if opts.tflab.loglevel > 0
        logmsg('Computing metrics on full time range of data using average Z and\n');
        logmsg('computing metrics on segments using segment Zs.\n');
    end
    S = tflab_metrics(S);
else

    logmsg('Computing Z using stack regression.\n');        

    S = tflab_preprocess(B,E,opts,'both',1,1);
    S = stackregression_(S,opts);

    logmsg('Computing metrics on full time series and segments using stack regression Z.\n');
    S = tflab_metrics(S);

    %logmsg('Computing metrics on segments using computed Z.\n');    
    %S = tflab_metrics(S,1);
end

end % tflab()

function S = stackregression_(S,opts)
    if opts.tflab.loglevel > 0
        logmsg('Computing Z using stack regression.\n');
    end
    DFT = dftcombine(S.Segment.DFT);
    [S.Z,S.fe,S.dZ,S.Regression] = tflab_miso(DFT,opts);
end

function S = stackaverage_(S,opts)

    logmsg('Computing Z by stack averaging segment Zs.\n');

    Sc = S.Segment;
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
        [Sc{s}.Z,Sc{s}.fe,Sc{s}.dZ,Sc{s}.Regression] = tflab_miso(Sc{s}.DFT,opts);

        if 0
            if opts.tflab.loglevel > 0
                logmsg('Computing metrics on segment using segment Z.\n');
            end
            Sc{s} = tflab_metrics(Sc{s});
            if opts.tflab.loglevel > 0
                logmsg('Computed metrics on segment using segment Z.\n');
            end
            if opts.tflab.loglevel > 0 && ~isempty(opts.fd.stack.average.function)
                logmsg('Computed Z and metrics for segment %d of %d; PE/CC/MSE %.2f/%.2f/%.3f\n',...
                       s,length(Sc),Sc{s}.Metrics.PE,Sc{s}.Metrics.CC,Sc{s}.Metrics.MSE);
            end
        end
    end

    S.Segment = combineStructs(Sc,3);
    S.Z = mean(S.Segment.Z,3);
    S.dZ = mean(S.Segment.dZ,3);
    S.fe = Sc{1}.fe;    

end

function S = intervals_(B, E, opts)
    
    if isempty(opts.fd.stack.average.function)
        S = tflab_preprocess(B,E,opts,'both',1,1);
        S = stackregression_(S,opts);
    else
        S = tflab_preprocess(B,E,opts,'both',1,0);
        S = stackaverage_(S,opts);
    end
    S = tflab_metrics(S);
end