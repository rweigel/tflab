function S = tflab_preprocess(In,Out,opts,domain,includeSegs,combineSegs)
%
%  S = TFLAB_PREPROCESS(S)
%
%  S = TFLAB_PREPROCESS(In, Out, opts)
%  S = TFLAB_PREPROCESS(In, Out, opts, domain)
%  S = TFLAB_PREPROCESS(In, Out, opts, domain, includeSegs)
%  S = TFLAB_PREPROCESS(In, Out, opts, domain, 1, combineSegs)

% _argcheck(varargin{:)

if isstruct(In)
    S = In;
    opts = S.Options;
    domain = 'both';
    combineSegs = 1;
    includeSegs = 1;
    if isfield(S,'Segment')
        % If running preprocess on saved state, saved state
        % may have results of regression. S.Segment will
        % be over-written and then fields of Segmento are recovered.
        Segmento = S.Segment;
        S = rmfield(S,'Segment');
    end
else
    if nargin < 3
        error('If In is a time series, at least 3 arguments are required.');
    end
    if nargin < 4
        domain = 'both';
    end
    if nargin < 5
        includeSegs = 0;
    end
    if nargin < 6
        combineSegs = 0;
    end
    if nargin == 6 && includeSegs ~= 1
        error('If combineSegs is given, includeSegs must be 1.');
    end
    if iscell(In)
        S = intervals_(In,Out,opts,domain,includeSegs,combineSegs);
        return;
    end
    S = struct('In',In,'Out',Out,'Options',opts);
end

if combineSegs && ~includeSegs
    % TODO: Fix awkward API
    error('If combineSegs = 1, includeSegs must be 1');
end

loglevel = opts.tflab.loglevel;

% Process full time series
S = tflab_tdpreprocess(S,loglevel);
S = tflab_fdpreprocess(S);

if ~isfield(S.Options,'td')
    return
end

Tw = opts.td.window.width;
Ts = opts.td.window.shift;

emsg = 'If one of window.width or window.shift is NaN, both must be NaN';
assert(any(isnan([Tw,Ts])) == all(isnan([Tw,Ts])),emsg);

if (isnan(Tw) && isnan(Ts)) || ~includeSegs %|| Tw == Ts
    % No segments to create or no segments desired.
    return
end

[a,b] = segmentidxs_(size(S.In,1),Tw,Ts);
for s = 1:length(a)

    Iseg = a(s):b(s);
    if loglevel > 0
        logmsg('Preprocessing segment %d of %d\n',s,length(a));
        logmsg('Segment time index range = [%d:%d]\n',Iseg(1),Iseg(end));
    end
    
    S.Segment{s} = struct();
    S.Segment{s}.Options = opts;
    S.Segment{s}.IndexRange = [a(s),b(s)]';

    S.Segment{s}.In  = S.In(Iseg,:);
    S.Segment{s}.Out = S.Out(Iseg,:);

    S.Segment{s} = tflab_tdpreprocess(S.Segment{s},loglevel);
    S.Segment{s} = tflab_fdpreprocess(S.Segment{s});
   
end

if combineSegs
    logmsg('Combining segments\n')
    S.Segment = combineStructs(S.Segment,3);
    if exist('Segmento','var')
        fns = fieldnames(Segmento);
        for i = 1:length(fns)
            % If running preprocess on saved state, saved state
            % may have results of regression. Here we recover them.
            S.Segment.(fns{i}) = Segmento.(fns{i});
        end
    end
    S.Segment = rmfield(S.Segment,'Options');
end


end % function

function [a,b] = segmentidxs_(N,Tw,Ts)
    
    if isnan(Tw) || isnan(Ts)
        a = 1;
        b = N;
        return
    end
    
    assert(Tw > 1, 'opts.td.window.width must be greater than 1');
    assert(Ts > 1, 'opts.td.window.shift must be greater than 1');

    assert(Tw <= N, 'opts.td.window.width must be <= size(In,1)');
    assert(mod(Tw,1) == 0, 'opts.td.window.width must be an integer');
    assert(mod(Ts,1) == 0, 'opts.td.window.shift must be an integer');        

    a = 1:Ts:N; % Start indices
    b = a + Tw - 1;      % Stop indices

    if b(end) > N
        % If last window not full, omit it.
        a = a(1:end-1);
        b = b(1:end-1);
    end

    if b(end) ~= N
        warning(['opts.td.window.width = %d,'...
                 'opts.td.window.shift = %d, size(In,1) = %d.'...
                 '\n\t Last %d point(s) will not be used.\n'],...
                 Tw, Ts, N, N-b(end));
    end
end

function S = intervals_(In, Out, opts, domain, includeSegs, combineSegs)    

    % Each cell contians an interval and an arbitrary gap in time is
    % assumed between the end of one cell and the start of the next
    % cell.
    assert(all(size(In) == size(Out)),'Required: size(B) == size(E)');
    assert(isvector(In),'Required: B must be vector cell array');
    assert(isvector(Out),'Required: E must be vector cell array');

    S.Options = opts;
    
    sOut = size(Out{1},2);
    sIn = size(In{1},2);

    % Check that number of columns is same for all intervals.
    for i = 2:length(In)
        assert(size(In{i},2) == sIn,...
            'Number of columns in In{%d} must be the same as Out{1}.',i);
        assert(size(Out{i},2) == sOut,...
            'Number of columns in In{%d} must be the same as Out{1}.',i);
        Nr(i) = size(In,1);
    end
    
    if isnan(opts.td.window.width) && length(unique(Nr)) ~= 1
        warning(['opts.td.window.width not given and intervals do not', ...
                'all have the same length.\n', ...
                'Using shortest interval length (%d) for width and shift'],...
                min(Nr));
        opts.td.window.width = min(Nr);
        opts.td.window.shift = min(Nr);
    end
    
    for i = 1:length(In) % Number of intervals
        Si{i} = tflab_preprocess(In{i},Out{i},opts,domain,includeSegs,combineSegs);        
        Si{i} = rmfield(Si{i},'Options');
    end
    p = 1;
    for i = 1:length(Si)
        fns = fieldnames(Si{i}); 
        for f = 1:length(fns)
            if ~any(strcmp(fns{f},{'Options','Segment'}))
                S.(fns{f}){i} = Si{i}.(fns{f});
            end
        end
        if combineSegs
            Segment{p} = Si{i}.Segment;
            p = p+1;            
        else
            for s = 1:length(Si{i}.Segment)
                Segment{p} = Si{i}.Segment{s};
                p = p+1;
            end
        end
    end
    for i = 1:length(Segment)
        if isfield(Segment{i},'Options')
            Segment{i} = rmfield(Segment{i},'Options');
        end
    end
    if combineSegs
        S.Segment = combineStructs(Segment,3);
    else
        S.Segment = Segment;        
    end

end