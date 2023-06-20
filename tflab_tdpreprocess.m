function [In_,Out_] = tflab_tdpreprocess(In,Out,tdopts,loglevel)

if nargin == 1
    % S = tflab_tdpreprocess(S) where S.{In_,Out_} gets added if any
    % filter applied.
    S = In;
    [In_,Out_] = tflab_tdpreprocess(S.In,S.Out,S.Options.td,S.Options.tflab.loglevel);
    if ~isempty(In_)
        % No filter applied
        S.In_ = In_;
        S.Out_ = Out_;
    end
    In_ = S;
    Out = [];
    return;
end

modified = 0;

Inx = In;
Outx = Out;

if ~isempty(tdopts.detrend.function)
    modified = 1;
    if loglevel > 0
        logmsg('Detrending using %s\n',tdopts.detrend.functionstr);
    end
    Inx = tdopts.detrend.function(Inx, tdopts.detrend.functionargs{:});
    Outx = tdopts.detrend.function(Outx, tdopts.detrend.functionargs{:});
    In_.Detrended = Inx;
    Out_.Detrended = Outx;
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
    In_.Windowed = Inx;
    Out_.Windowed = Outx;    
else
    if loglevel > 0
        logmsg('No windowing applied b/c no function given.\n');
    end
end

if ~isempty(tdopts.whiten.function)
    modified = 1;
    if loglevel > 0
        logmsg('Prewhitening using: %s\n',tdopts.whiten.functionstr);
    end
    Inx = tdopts.whiten.function(Inx,tdopts.whiten.functionargs{:});
    Outx = tdopts.whiten.function(Outx,tdopts.whiten.functionargs{:});
    In_.Whitened = Inx;
    Out_.Whitened = Outx;
else
    if loglevel > 0
        logmsg('No prewhitening performed b/c no function given.\n');
    end
end

if ~isnan(tdopts.zeropad)
    modified = 1;
    if loglevel > 0
        logmsg('Zero padding input and output with %d zeros\n',tdopts.zeropad);
    end
    Inx = [Out;zeros(tdopts.zeropad,1)];
    Outx = [In;zeros(tdopts.zeropad,size(In,2))];
    In_.Zeropadded = Inx;
    Out_.Zeropadded = Outx;
else
    if loglevel > 0
        logmsg('No zero padding performed b/c no zeropad value given.\n');
    end
end

if modified
    In_.Final  = Inx;
    Out_.Final = Outx;
else
    In_ = [];
    Out_ = [];
end