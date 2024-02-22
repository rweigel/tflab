function [Z,fe,dZ,Regression] = tflab_miso(DFT,opts)
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
    logmsg('Doing regression using %s\n',opts.fd.regression.functionstr);
end

boot_note = 1;
dZ = [];
for j = 1:length(fe)

    ftIn = DFTIn{j,1};
    ftOut = DFTOut{j,1};

    Z(j,:) = nan(1,size(ftIn,2));
    dZ(j,1) = nan;
    Regression.Residuals{j,1} = nan(size(ftOut,1),1);

    if opts.fd.window.loglevel > 0
        msg = 'Band with center of fe = %.8f has %d points; fl = %.8f fh = %.8f\n';
        logmsg(msg,fe(j),length(f),min(f),max(f));
    end

    if 0 && size(ftIn,2) == 1 && length(f) == 1
        % One input component
        z = ftOut./ftIn;
        if isinf(z)
            z = nan;
        end
        Z(j,1) = z;
        dZ(j,1) = 0;
        continue;
    end
    
    if length(f{j}) < 2*size(ftIn,2)
        msg = '!!! System is underdetermined for fe = %f. Setting Z equal to NaN(s).\n';
        logmsg(msg,fe(j));
        continue;
    end

    if length(f{j}) == 2*size(ftIn,2)
        msg = '!!! System is exactly determined for fe = %f. Setting Z equal to NaN(s).\n';
        logmsg(msg,fe(j));
        continue;
    end

    lastwarn('');

    regressargs = opts.fd.regression.functionargs;
    regressfunc = opts.fd.regression.function;
    [Z(j,:),dz,Info] = regressfunc(ftOut,ftIn,regressargs{:});

    if ~isempty(dz)
        dZ(j,1) = dz;
    end

    if ~isempty(lastwarn)
        msg = 'Above warning is for eval. freq. #%d; fe = %f; Te = %f\n';
        logmsg(msg,j,fe(j),1/fe(j));
        logmsg('ftE =');
        ftOut
        logmsg('ftB =');
        ftIn
    end

    if any(isinf(Z(j,:)))
        msg = '!!! Z has Infs for fe = %f. Setting all element of Z to NaN(s).\n';
        logmsg(msg,fe(j));
        Z(j,:) = nan;
        dZ(j,1) = nan;
        continue;
    end

    if isfield(Info,'Residuals')
        % Residuals are also computed in tflab_metrics because we remove
        % Regression.Residuals when saving file to reduce file size. The
        % following is not needed but is kept because if tflab_metrics
        % finds this, it will check that its calculation matches.
        Regression.Residuals{j,1} = Info.Residuals;
    end
    if isfield(Info,'ZCL95l')
        Regression.ZCL95l(j,:) = Info.ZCL95l;
    end
    if isfield(Info,'ZCL95u')
        Regression.ZCL95u(j,:) = Info.ZCL95u;
    end
    if isfield(Info,'dZCL95l')
        Regression.dZCL95l(j,:) = Info.dZCL95l;
    end
    if isfield(Info,'dZCL95u')
        Regression.dZCL95u(j,:) = Info.dZCL95u;
    end

    n = size(ftOut,1);
    if opts.fd.zerrorbars_ && n > 10
        Nb = 100; % Number of bootstrap samples
        fract = 1; % Fraction to sample; Efron's original bootstrap method uses 1;
                   % If f != 1, called m of n bootstrap; see 10.1002/9781118445112.stat08002
        if boot_note == 1
            logmsg(sprintf('Computing confidence limits using %d bootstrap samples and m/n = %.2f\n',Nb,fract));
            boot_note = 0;
        end
        for b = 1:Nb
            I = randsample(n,round(fract*n),1); % Resample with replacement
            Zb(b,:) = regressfunc(ftOut(I,:),ftIn(I,:),regressargs{:});
        end
        nl = round((0.05/2)*Nb);
        nh = round((1-0.05/2)*Nb);

        for c = 1:size(Z,2)
            Zb(:,c) = sort(abs(Zb(:,c)),1); % Sort rows
            l = abs(Zb(nl,c));    % Select the nth lowest
            u = abs(Zb(nh,c));    % Select the nth highest
            Regression.Bootstrap.ZVAR(j,c) = var(abs(Zb(:,c)),0,1);
            Regression.Bootstrap.ZCL95l(j,c) = l;
            Regression.Bootstrap.ZCL95u(j,c) = u;
        end
    end

end

if opts.fd.regression.loglevel > 0
    logmsg(['Finished freq band and regression '...
            'calculations for %d eval. freqs.\n'],...
             length(fe)-1);
end

if all(isnan(Z(:)))
    error('All Z values are NaN');
end

