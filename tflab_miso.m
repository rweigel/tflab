function [Z,fe,dZ,Regression] = tflab_miso(DFT,opts,error_estimates)
%TFLAB_MISO Frequency domain MISO transfer function estimate
%
%  S = TFLAB_MISO(DFT) returns a structure with an estimate of the
%  transfer function Z in the expression
%
%    E(f) = Zx(f)Bx(f) + Zy(f)By(f) + ...
%
%  given time series for E(t), Bx(t), By(t), ...
%
%  Estimates are made for the complex-valued transfer function Z at a set of
%  evaluation frequencies, fe, using regression with a set of frequencies in
%  a band around each fe on the model equation above.
%
%  See also TFLAB.

Regression = struct();

if isfield(DFT,'In_')
    msg = 'Using DFTs from filtered In and Out (DFT.In_.Final and DFT.Out_.Final)';
    logmsg(msg);
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
    msg = 'Starting freq band and regression calcs for %d frequencies.\n';
    logmsg(msg,length(fe));
    logmsg('Doing regression using %s\n',opts.fd.regression.functionstr);
end

boot_note = 1;
dZ = [];
for j = 1:length(fe)

    ftIn = DFTIn{j,1};
    ftOut = DFTOut{j,1};

    Z(j,:) = (1+1j)*nan(1,size(ftIn,2));
    dZ(j,1) = (1+1j)*nan;
    Residuals{j,1} = (1+1j)*nan(size(ftOut,1),1);
    Parametric.ZCL95l(j,:)  = (1+1j)*nan*ones(1,size(Z,2));
    Parametric.ZCL95u(j,:)  = (1+1j)*nan*ones(1,size(Z,2));
    Parametric.dZCL95l(j,1) = (1+1j)*nan;
    Parametric.dZCL95u(j,1) = (1+1j)*nan;

    if opts.fd.window.loglevel > 0
        msg = 'Band with center of fe = %.8f has %d points; fl = %.8f fh = %.8f\n';
        logmsg(msg,fe(j),length(f),min(f),max(f));
    end

    if 0 && size(ftIn,2) == 1 && length(f) == 1
        % One input component
        z = ftOut./ftIn;
        if isinf(z)
            z = nan*(1+1j);
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
        Z(j,:) = (1+1j)*nan;
        dZ(j,1) = (1+1j)*nan;
        continue
    end

    if isfield(Info,'Residuals')
        % Residuals are also computed in tflab_metrics because we remove
        % Regression.Residuals when saving file to reduce file size. The
        % following is not needed but is kept because if tflab_metrics
        % finds this, it will check that its calculation matches.
        Residuals{j,1} = Info.Residuals;
    end
    if error_estimates
        % Don't keep error estimats if transfer function is computed based on
        % segment stack averages.
        if isfield(Info,'ZCL95l')
            Parametric.ZCL95l(j,:) = Info.ZCL95l;
        end
        if isfield(Info,'ZCL95u')
            Parametric.ZCL95u(j,:) = Info.ZCL95u;
        end
        if isfield(Info,'dZCL95l')
            Parametric.dZCL95l(j,:) = Info.dZCL95l;
        end
        if isfield(Info,'dZCL95u')
            Parametric.dZCL95u(j,:) = Info.dZCL95u;
        end
    end

    if ~isfield_(opts,'bootstrap.fd') && j == 1
        logmsg('Not computing bootstrap error estimates b/c/ boostrap.fd options not given.\n');
        continue
    end
    if ~error_estimates && j == 1
        logmsg('Bootstrap error estimates not requested by calling function.\n');
        continue
    end
    n = size(ftOut,1);
    if n >= opts.fd.bootstrap.nmin
        % Bootstrap confidence limits if requested, nmin samples or more, and
        % transfer function error estimates are not computed a different way.
        Nb = opts.fd.bootstrap.N;
        fraction = opts.fd.bootstrap.fraction;
        m = round(fraction*n);
        if boot_note == 1
            msg = 'Computing confidence limits using %d bootstrap samples and m/n = %.2f\n';
            logmsg(msg,Nb,fraction);
            boot_note = 0;
        end
        for b = 1:Nb
            I = randsample(n,m,1); % Resample with replacement
            Zb(b,:) = regressfunc(ftOut(I,:),ftIn(I,:),regressargs{:});
        end
        Bootstrap(j) = error_estimates_bootstrap(fe(j),Zb,Z(j,:));
    end
end

Regression.Residuals = Residuals;
Regression.ErrorEstimates.Parametric = Parametric;
if exist('Bootstrap','var')
    Regression.ErrorEstimates.Bootstrap = combineBootstrap(Bootstrap);
end

if opts.fd.regression.loglevel > 0
    msg = 'Finished freq band and regression calculations for %d eval. freqs.\n';
    logmsg(msg,length(fe)-1);
end

if all(isnan(Z(:)))
    error('All Z values are NaN');
end

end % function tflab_miso()
