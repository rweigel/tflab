function S = tflab_preprocess(S);

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

B = S.In;
E = S.Out;

if ~isempty(opts.td.detrend.function)
    if opts.tflab.loglevel > 0
        logmsg('Detrending input and output using %s\n',...
                opts.td.detrend.functionstr);
    end
    B = opts.td.detrend.function(B, opts.td.detrend.functionargs{:});
    E = opts.td.detrend.function(E, opts.td.detrend.functionargs{:});
    S.Detrend.In = B;
    S.Detrend.Out = E;
else
    if opts.tflab.loglevel > 0
        logmsg('No detrending b/c no function given.\n');
    end
end

% TODO?: Allow TD window and prewhiten to not be same for input and output
% and then compute corrected Z.
if ~isempty(opts.td.window.function)
    if opts.tflab.loglevel > 0
        logmsg('Windowing input and output using %s\n',...
                opts.td.window.functionstr);
    end
    [B,~] = opts.td.window.function(B,opts.td.window.functionargs{:});
    [E,W] = opts.td.window.function(E,opts.td.window.functionargs{:});
    S.Window.Weights = W;
    S.Window.In = B;
    S.Window.Out = E;
else
    if opts.tflab.loglevel > 0
        logmsg('No time domain window applied b/c no function given.\n');
    end
end

if ~isempty(opts.td.prewhiten.function)
    if opts.tflab.loglevel > 0
        logmsg('Prewhitening input and output using: %s\n',...
                opts.td.prewhiten.functionstr);
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
        logmsg(['No time domain prewhitening performed '...
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
    logmsg( 'Computing DFT of preprocessed input and output(s).\n');
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

winfn = opts.fd.window.function;
if opts.fd.window.loglevel && strcmp(opts.fd.program.name,'tflab')
    logmsg( 'Using FD window function %s\n',func2str(winfn));
end

for j = 1:length(Ic)

    %logmsg('Computing DFT for freq. %d of %d\n',j, length(fe));

    W = winfn(2*Ne(j)+1);
    W = W/sum(W);
    W  = sqrt(W);

    r = Ic(j)-Ne(j):Ic(j)+Ne(j); % Index range

    S.DFT.Out{j,1} = ftE(r,1);
    S.DFT.In{j,1} = ftB(r,:);
    S.DFT.f{j,1} = f(r);
    S.DFT.Weights{j,1} = W;
end

%[S.Metrics.PSD.Raw.In,S.Metrics.DFT.Raw.In,S.Metrics.DFT.Raw.fe] = psd(B);
%[S.Metrics.PSD.Raw.Out,S.Metrics.DFT.Raw.Out] = psd(E);
