function S = tflab_miso(B,E,t,opts)
%TFLAB_MISO Frequency domain MISO transfer function estimate
%
%  S = TFLAB_MISO(B,E) returns a structure with an estimate of the
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
%  opts.td.window.width and opts.td.window.shift are ignored.
%
%  See also TFLAB.

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
        logmsg('Prewhitening input and output using: %s\n',opts.td.prewhiten.functionstr);
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
                %logmsg(['System is underdetermined for fe = %f. ',...
                %        'Setting Z equal to zero(s) for this frequency.'],fe(j));
            else
                Z(j,:) = nan(1,size(ftB,2));
                S.Regression.Weights{j,1} = nan*W;
                S.Regression.Residuals{j,1} = nan*W;
                logmsg(['!!! System is underdetermined for fe = %f. ',...
                        'Setting Z equal to NaN(s) for this frequency.'],fe(j));
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
                logmsg('!!! System if underdetermined for fe = %f. Setting Z equal to zero(s) for this frequency.',fe(j));
                continue;
            end
        end

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