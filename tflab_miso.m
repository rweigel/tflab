function [Z,fe,Regression] = tflab_miso(DFT,opts)
%TFLAB_MISO Frequency domain MISO transfer function estimate
%
%  S = TFLAB_MISO(DFT) returns a structure with an estimate of the
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

Regression = struct();

%function [Z,fe,Regression] = tflab_miso(DFTIn, DFTOut, f, fe, opts)

if isfield(DFT,'In_')
    logmsg('Using DFTs from filtered In and Out (DFT.In_.Final and DFT.Out_.Final)');
    DFTIn = DFT.In_.Final;
    DFTOut = DFT.Out_.Final;
    f = DFT.f_;
    fe = DFT.fe_;
else
    DFTIn = DFT.In;
    DFTOut = DFT.Out;
    f = DFT.f;
    fe = DFT.fe;
end

if opts.tflab.loglevel > 0
    logmsg(['Starting freq band and regression '...
            'calcs for %d frequencies.\n'],length(fe));
    logmsg(['Using %s() with additional arguments given in\n'...
            'opts.fd.regression.functionargs\n'],...
            func2str(opts.fd.regression.function));
end

for j = 1:length(fe)

    ftIn = DFTIn{j,1};
    ftOut = DFTOut{j,1};

    if opts.fd.window.loglevel > 0
        logmsg(['Band with center of fe = %.8f has %d '...
                'points; fl = %.8f fh = %.8f\n'],...
                 fe(j),length(f),min(f),max(f));
    end

    if size(ftIn,2) == 1 && length(f) == 1
        z = ftOut./ftIn;
        if isinf(z)
            z = nan;
        end
        keyboard
        Z(j,1) = z;
        Regression.Weights{j,1} = nan(length(f),1);
        Regression.Residuals{j,1} = nan(length(f),1);
        continue;
    end
    
    if length(f) < size(ftIn,2)
        Z(j,:) = nan(1,size(ftIn,2));
        Regression.Weights{j,1} = nan(length(f),1);
        Regression.Residuals{j,1} = nan(length(f),1);
        logmsg(['!!! System is underdetermined for fe = %f. ',...
                'Setting Z equal to NaN(s) for this frequency.'],fe(j));
        continue;
    end

    lastwarn('');
    regressargs = opts.fd.regression.functionargs;
    regressfunc = opts.fd.regression.function;
    
    [Z(j,:),Residuals,Weights] = ...
               regressfunc(ftOut,ftIn,regressargs{:});
    
    n = size(ftOut,1);
    if  n > 10
        Nb = 1000;
        for b = 1:Nb
            I = randsample(n,round(0.63*n),1);
            Zb(b,:) = regressfunc(ftOut(I,:),ftIn(I,:),regressargs{:});
        end
        nl = round((1-0.95)*Nb);
        nh = round(0.95*Nb);
        Zb = sort(abs(Zb),1); % Sort
        l = Zb(nl,:);    % Select the nth lowest
        u = Zb(nh,:);    % Select the nth highest
        Regression.ZVAR(j,:) = var(Zb,0,1);
        Regression.ZCL(j,1:2:2*length(l)) = l;
        Regression.ZCL(j,2:2:2*length(u)+1) = u;
    end

    if ~isempty(lastwarn)
        logmsg('Above is for eval. freq. #%d; fe = %f; Te = %f\n', ...
            j,fe(j),1/fe(j));
        logmsg(sprintf('ftE = \n'));
        logmsg(sprintf('   %.16f\n',ftOut));
        logmsg(sprintf('ftB = \n'));
        logmsg(sprintf('   %.16f\n',ftIn));
    end

    if any(isinf(Z(j,:)))
        logmsg(['!!! Z has Infs for fe = %f. ',...
                'Setting Z equal to NaN(s) for this frequency.'],fe(j));
        Z(j,:) = nan;
        Regression.Weights{j,1} = nan(length(f),1);
        Regression.Residuals{j,1} = nan(length(f),1);
        continue;
    end
           
    Regression.Weights{j,1} = Weights;
    Regression.Residuals{j,1} = Residuals;

end

if opts.fd.regression.loglevel > 0
    logmsg(['Finished freq band and regression '...
            'calculations for %d eval. freqs.\n'],...
             length(fe)-1);
end

if all(isnan(Z(:)))
    error('All Z values are NaN');
end

