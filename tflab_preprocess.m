function S = tflab_preprocess(In,Out,opts)

if nargin == 1
    % tflab_preprocess(S)
    S = In;
    optso = S.Options;
    optsx = optso;
    optsx.td.window.width = NaN;
    optsx.td.window.shift = NaN;
    % Do calcs for full time series
    S = tflab_preprocess(S.In, S.Out, optsx);
    if ~isnan(optso.td.window.width)
        % Do calcs for segments
        Segments = tflab_preprocess(S.In, S.Out, optso);
        S.Segment = combineStructs(Segments,3);
        S.Segment = rmfield(S.Segment,'Options');
    end
    return
end

Tw = opts.td.window.width;
Ts = opts.td.window.shift;
assert(any(isnan([Tw,Ts])) == all(isnan([Tw,Ts])),...
        'If one of window.width or window.shift is NaN, both must be NaN');
    
[a,b] = segmentidxs_(size(In,1),Tw,Ts);

for s = 1:length(a)

    Iseg = a(s):b(s);
    if opts.tflab.loglevel > 0
        logmsg('Preprocessing segment %d of %d\n',s,length(a));
        logmsg('Segment time index range = [%d:%d]\n',Iseg(1),Iseg(end));
    end
    
    S{s} = struct();
    S{s}.Options = opts;
    S{s}.IndexRange = [a(s),b(s)]';

    S{s}.In  = In(Iseg,:);
    S{s}.Out = Out(Iseg);
    
    S{s} = tflab_tdpreprocess(S{s});
    S{s} = tflab_fdpreprocess(S{s});

end
if length(S) == 1
    S = S{1};
    S = rmfield(S,'IndexRange');
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