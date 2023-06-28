function S = tflab_tdpreprocess(S,loglevel)

if nargin < 2
    loglevel = 0;
end

tdopts = S.Options.td;

modified = 0;

Inx = S.In;
Outx = S.Out;

if ~isempty(tdopts.detrend.function)
    modified = 1;
    if loglevel > 0
        logmsg('Detrending using %s\n',tdopts.detrend.functionstr);
    end
    Inx = tdopts.detrend.function(Inx, tdopts.detrend.functionargs{:});
    Outx = tdopts.detrend.function(Outx, tdopts.detrend.functionargs{:});
    S.In_.Detrended = Inx;
    S.Out_.Detrended = Outx;
else
    if loglevel > 0
        logmsg('No detrending b/c no function given.\n');
    end
end

% TODO?: Allow TD window and prewhiten to not be same for input and output
% and then compute corrected Z.
if ~isempty(tdopts.window.function)
    modified = 1;
    if loglevel > 0
        logmsg('Windowing using %s\n',tdopts.window.functionstr);
    end
    Inx = tdopts.window.function(Inx,tdopts.window.functionargs{:});
    Outx = tdopts.window.function(Outx,tdopts.window.functionargs{:});
    S.In_.Windowed = Inx;
    S.Out_.Windowed = Outx;    
else
    if loglevel > 0
        logmsg('No windowing applied b/c no window function given.\n');
    end
end

if ~isempty(tdopts.whiten.function)
    modified = 1;
    if loglevel > 0
        logmsg('Whitening using: %s\n',tdopts.whiten.functionstr);
    end
    Inx = tdopts.whiten.function(Inx,tdopts.whiten.functionargs{:});
    Outx = tdopts.whiten.function(Outx,tdopts.whiten.functionargs{:});
    S.In_.Whitened = Inx;
    S.Out_.Whitened = Outx;
else
    if loglevel > 0
        logmsg('Whitening performed b/c no whiten function given.\n');
    end
end

if ~isnan(tdopts.zeropad)
    modified = 1;
    if loglevel > 0
        logmsg('Zero padding input and output with %d zeros\n',tdopts.zeropad);
    end
    Inx = [Out;zeros(tdopts.zeropad,1)];
    Outx = [In;zeros(tdopts.zeropad,size(In,2))];
    S.In_.Zeropadded = Inx;
    S.Out_.Zeropadded = Outx;
else
    if loglevel > 0
        logmsg('No zero padding performed b/c no zeropad value given.\n');
    end
end

if modified
    S.In_.Final  = Inx;
    S.Out_.Final = Outx;
end