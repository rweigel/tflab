function [Z,fe,Regression] = tflab_miso(DFT,opts,error_estimates)
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
for j = 1:length(fe)

    ftOut = DFTOut{j,1};
    ftIn = DFTIn{j,1};
    if opts.fd.regression.const_term
        ftIn = [ftIn,ones(size(ftIn,1),1)];
    end

    Z(j,:) = (1+1j)*nan(1,size(ftIn,2));
    Residuals{j,1} = (1+1j)*nan(size(ftOut,1),1);
    Parametric.ZCL95l(j,:)  = (1+1j)*nan*ones(1,size(Z,2));
    Parametric.ZCL95u(j,:)  = (1+1j)*nan*ones(1,size(Z,2));

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
        continue;
    end

    if length(f{j}) < 2*(size(ftIn,2))
        msg = '!!! System is underdetermined for fe = %.1e. Setting Z equal to NaN(s).\n';
        logmsg(msg,fe(j));
        continue;
    end

    if length(f{j}) == 2*(size(ftIn,2))
        msg = '!!! System is exactly determined for fe = %.1e. Setting Z equal to NaN(s).\n';
        logmsg(msg,fe(j));
        continue;
    end

    [Z(j,:),Info] = callregress_(ftOut,ftIn,j,fe(j),opts);

    if any(isinf(Z(j,:)))
        continue;
    end

    if isfield(Info,'Residuals')
        % Residuals are also computed in tflab_metrics because we remove
        % Regression.Residuals when saving file to reduce file size. The
        % following is not needed but is kept because if tflab_metrics
        % finds this, it will check that its calculation matches.
        Residuals{j,1} = Info.Residuals;
    end

    if error_estimates
        % error_estimates = 0 if transfer function is computed based on
        % segment stack averages. In that case, we don't keep error estimates
        % for the Z computed for each segment.
        if isfield(Info,'ZCL95l')
            Parametric.ZCL95l(j,:) = Info.ZCL95l;
        end
        if isfield(Info,'ZCL95u')
            Parametric.ZCL95u(j,:) = Info.ZCL95u;
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
    if ~isnan(opts.fd.bootstrap.N) && n >= opts.fd.bootstrap.nmin
        % Bootstrap confidence limits if requested, nmin samples or more, and
        % transfer function error estimates are not computed a different way.
        Nb = opts.fd.bootstrap.N;
        fraction = opts.fd.bootstrap.fraction;
        m = round(fraction*n);
        if boot_note == 1
            msg = 'Computing Z confidence limits using %d bootstrap samples and m/n = %.2f\n';
            logmsg(msg,Nb,fraction);
            boot_note = 0;
        end
        Zb = [];
        parfor b = 1:Nb
            I = randsample(n,m,1); % Resample with replacement
            if length(unique(I)) < 2*(size(ftIn,2) + 1)
                msg = '!!! System is underdetermined for bootstrap resample %d for fe = %.1e.\n';
                logmsg(msg,b,fe(j));
                msg = '!!! Skipping this resampling. Consider increasing bootstrap.nmin to be >> 2*size(In,2) + 1.\n';
                logmsg(msg);
                continue;
            end
            Zb = [Zb; callregress_(ftOut(I,:),ftIn(I,:),j,fe(j),opts)];
        end
        if size(Zb,1) >= opts.fd.bootstrap.nmin
            Bootstrap(j) = error_estimates_bootstrap(fe(j),Zb,Z(j,:));
        else
            msg1 = 'cannot compute boostrap error estimates b/c ';
            msg2 = sprintf('# of bootstrap samples (%d) < opts.fd.bootstrap.nmin (%d)',...
                           size(Zb,1),opts.fd.bootstrap.nmin);
            logmsg("!!! For f = %.1e, %s%s due to sample skipping.\n", fe(j), msg1, msg2);
        end
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

function [Z,Info] = callregress_(ftOut,ftIn,j,fe,opts)

    regressargs = opts.fd.regression.functionargs;
    regressfunc = opts.fd.regression.function;

    warning('off','stats:regress:NoConst'); % Suppress display of warning.
    lastwarn('');

    [Z,Info] = regressfunc(ftOut,ftIn,regressargs{:});

    % Ideally we would use the following to check if the warning was
    % for stats:regress:NoConst, but it doesn't work. So we use a string
    % match on lastwarn to determine if the warning was for NoConst.
    %  warning('on','stats:regress:NoConst');
    %  w = warning('query','last');
    %  strcmp(w.identifier,'stats:regress:NoConst')

    msg = 'R-square and the F statistic are not well-defined';
    if ~isempty(lastwarn) && ~startsWith(lastwarn, msg)
        msg = 'Above warning is for eval. freq. #%d; fe = %f; Te = %f\n';
        logmsg(msg,j,fe,1/fe);
        logmsg('ftE =');
        ftOut
        logmsg('ftB =');
        ftIn
    end

    if any(isinf(Z))
        msg = '!!! Z has Infs for fe = %f. Setting all elements of Z to NaN(s).\n';
        logmsg(msg,fe);
        Z = (1+1j)*nan*ones(1,size(Z,2));
    end
end % function callregress_()
