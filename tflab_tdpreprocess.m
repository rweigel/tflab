function S = tflab_tdpreprocess(S)

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

B = S.In;
E = S.Out;

if ~isempty(opts.td.detrend.function)
    if opts.tflab.loglevel > 0
        logmsg('Detrending using %s\n',opts.td.detrend.functionstr);
    end
    B = opts.td.detrend.function(B, opts.td.detrend.functionargs{:});
    E = opts.td.detrend.function(E, opts.td.detrend.functionargs{:});
    S.In_.Detrended = B;
    S.Out_.Detrended = E;
else
    if opts.tflab.loglevel > 0
        logmsg('No detrending b/c no function given.\n');
    end
end

% TODO?: Allow TD window and prewhiten to not be same for input and output
% and then compute corrected Z.
if ~isempty(opts.td.window.function)
    if opts.tflab.loglevel > 0
        logmsg('Windowing using %s\n',opts.td.window.functionstr);
    end
    B = opts.td.window.function(B,opts.td.window.functionargs{:});
    E = opts.td.window.function(E,opts.td.window.functionargs{:});
    S.In_.Windowed = B;
    S.Out_.Windowed = E;
else
    if opts.tflab.loglevel > 0
        logmsg('No windowing applied b/c no function given.\n');
    end
end

if ~isempty(opts.td.prewhiten.function)
    if opts.tflab.loglevel > 0
        logmsg('Prewhitening using: %s\n',opts.td.prewhiten.functionstr);
    end
    B = opts.td.prewhiten.function(B,opts.td.prewhiten.functionargs{:});
    E = opts.td.prewhiten.function(E,opts.td.prewhiten.functionargs{:});
    S.In_.Whitened = B;
    S.Out_.Whitened = E;
else
    if opts.td.prewhiten.loglevel
        logmsg('No prewhitening performed b/c no function given');
    end
end

if ~isnan(opts.td.zeropad)
    if opts.tflab.loglevel > 0
        logmsg('Zero padding input and output with %d zeros\n',opts.td.zeropad);
    end
    E = [E;zeros(opts.td.zeropad,1)];
    B = [B;zeros(opts.td.zeropad,size(B,2))];
    S.In_.Zeropadded = B;
    S.Out_.Zeropadded = E;
else
    if opts.td.prewhiten.loglevel
        logmsg('No zero padding performed b/c values given.');
    end
end

S.In_.Final  = B;
S.Out_.Final = E;

