function S = tflab_miso(S)
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

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

S.Regression = struct();

fe = S.fe;

if opts.tflab.loglevel > 0
    logmsg(['Starting freq band and regression '...
            'calcs for %d frequencies.\n'],length(fe));
    logmsg(['Using %s() with additional arguments given in\n'...
            'opts.fd.regression.functionargs\n'],...
            func2str(opts.fd.regression.function));
end

for j = 1:length(fe)

    ftB = S.DFT.In{j,1};
    ftE = S.DFT.Out{j,1};
    f = S.DFT.f{j,1};
    W = S.DFT.Weights{j,1};
    Wr = repmat(W,1,size(ftB,2));

    if opts.fd.window.loglevel > 0
        logmsg(['Band with center of fe = %.8f has %d '...
                'points; fl = %.8f fh = %.8f\n'],...
                 fe(j),length(f),f(1),f(end));
    end

    regressargs = opts.fd.regression.functionargs;
    regressfunc = opts.fd.regression.function;

    lastwarn('');

    if length(f) < size(ftB,2)
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

    if size(ftB,1) == 1 && size(ftB,2) == 1
        if ftB == 0 && ftE ~= 0 && fe(j) == 0.5
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
               regressfunc(W.*ftE,Wr.*ftB,regressargs{:});

    if any(isinf(Z(j,:)))
        logmsg(['!!! Z has Infs for fe = %f. ',...
                'Setting Z equal to NaN(s) for this frequency.'],fe(j));
        Z(j,:) = nan;
    end
    if ~isempty(lastwarn)
        logmsg('Above is for eval. freq. #%d; fe = %f; Te = %f\n', ...
            j,fe(j),1/fe(j));
        logmsg(sprintf('ftE = \n'));
        logmsg(sprintf('   %.16f\n',ftE));
        logmsg(sprintf('ftB = \n'));
        logmsg(sprintf('   %.16f\n',ftB));
    end

    S.Regression.Weights{j,1} = Weights;
    S.Regression.Residuals{j,1} = Residuals;

end

if opts.fd.regression.loglevel > 0
    logmsg(['Finished freq band and regression '...
            'calculations for %d eval. freqs.\n'],...
             length(Ic)-1);
end

if all(isnan(Z(:)))
    error('All Z values are NaN');
end

S.Z = Z;
S.Phi = atan2(imag(S.Z),real(S.Z));
