function ErrorEstimates = error_estimates(S,onsegments)

if onsegments
    Z = S.Segment.Z;
    fe = S.Segment.fe;
    t = tinv(0.975,size(Z,3)-1);
    for comp = 1:size(Z,2) % Second dimension is component
        Zc = squeeze(Z(:,comp,:)); % Columns of Zc are segment Zs
        STDr = std(real(Zc),1,2);
        STDi = std(imag(Zc),1,2);
        delta = t*(STDr + 1j*STDi)/sqrt(size(Z,3));
        Zave = mean(Zc,2);
        ZCL95l(:,comp) = Zave - delta;
        ZCL95u(:,comp) = Zave + delta;
    end
    ErrorEstimates.Parametric = error_estimates_derived(fe,Z,ZCL95l,ZCL95u);

    if isfield_(S,'Options.fd.bootstrap')
        if S.Options.fd.bootstrap.nmin < size(Z,3)
            logmsg('Not enough segments for bootstrap error estimates on segment Zs.\n')
            return
        end
        logmsg('Computing bootstrap error estimates on segment Zs.\n')
        for j = 1:size(S.Segment.Z,1)
            Zb = reshape(transpose(squeeze(S.Segment.Z(j,:,:))),size(Z,3),size(Z,2));
            Bootstrap(j) = error_estimates_bootstrap(fe(j),Zb,S.Z(j,:));
        end
        ErrorEstimates.Bootstrap = combineBootstrap(Bootstrap);
    end
    return
end

%logmsg('Computing additional error estimates.\n');

if isfield_(S,'Regression.ErrorEstimates.Bootstrap')
    % The Regression field is removed on save. So move into ErrorEstimates.
    ErrorEstimates.Bootstrap = S.Regression.ErrorEstimates.Bootstrap;
end

if isfield_(S,'Regression.ErrorEstimates.Parametric')
    logmsg('Computing derived parametric CLs based on ZCL95{u,l}.\n');

    % ZCL95{l,u} are (complex-valued) error estimates on Z, not |Z|.
    ZCL95l = S.Regression.ErrorEstimates.Parametric.ZCL95l;
    ZCL95u = S.Regression.ErrorEstimates.Parametric.ZCL95u;
    ErrorEstimates.Parametric = error_estimates_derived(S.fe,S.Z,ZCL95l,ZCL95u);
    return;
else
    if ~isfield(S,'ZVAR')
        logmsg('No ZVAR field at top-level and no ZCL95l in Regression.ErrorEstimates.Parametric.\n');
        logmsg('No error estimates available.\n');
        return
    end
    for comp = 1:size(S.Z,2)

        % This will occcur when S is based on content of EDI file, so no
        % TFLab regression information is available.

        logmsg('Found ZVAR field at top-level and no ZCL95l in Regression.ErrorEstimates.\n');
        logmsg('Computing ZMAGCL95{u,l} using it.\n');

        ZVAR(:,comp) = S.ZVAR(:,comp);

        sf = t*(1+1j)/sqrt(2);
        ZCL95l(:,comp) = S.Z(:,comp) - sf*sqrt(ZVAR(:,comp));
        ZCL95u(:,comp) = S.Z(:,comp) + sf*sqrt(ZVAR(:,comp));
    end
    ErrorEstimates.Parametric = error_estimates_derived(fe,Z,ZCL95l,ZCL95u);
end

end % function error_estimates()


% 1. Computing error bars from SPUD XML ZVAR.
%
% From https://ds.iris.edu/spud/resources/pdf/SPUD-XML-change-log.pdf
% (cached in runs/data/EarthScope/SPUD-XML-change-log.pdf)
% "To compute the standard error of a real or imaginary part, e.g.,
% for plotting, one would divide the variance by 2, then take square root
% (i.e. Std = sqrt(Var/2)). To plot the error bars for a real or imaginary part,
% one would then multiply that number by 2 ..."
%
% Using ZVAR from http://ds.iris.edu/spud/emtf/15014571, to compute delta_Z
% to use in RHOVAR calculation, I need to multiply delta_Z by sqrt(2), to
% get an approxmate visual match to the error bars for RHO shown (larger
% image at http://ds.iris.edu/spud/emtf/15014571). The documentation does
% not indicate if the error bars are 1 or 2 Standard Errors. Also, it may
% be that the plots have not been updated since 2018 when the reprocessing
% of the files was performed (see
% https://ds.iris.edu/spud/resources/pdf/SPUD-XML-change-log.pdf).
%
% TODO: Find out what is plotted at http://ds.iris.edu/spud/emtf/15014571
% so that I can state the error bar length in SEs.
%
% See also
%   https://ds.iris.edu/files/products/emtf/definitions/
%   http://ds.iris.edu/files/products/emtf/definitions/variance.html
%   https://library.seg.org/doi/epdf/10.1190/geo2018-0679.1
%
%
% 2. From EDI documentation (https://www.mtnet.info/docs/seg_mt_emap_1987.pdf),
%
% In order to provide rigorous error estimates for the surface impedance
% tensor elements, one must calculate the variance for the real part, the
% variance for the imaginary part, and a covariance. This is done in the
% usual manner for a finite number of estimates of the real and imaginary
% parts, which are clearly independent functions. This defines an error
% "ellipse" about the impedance tensor element in the complex plane. The
% >Z**R.VAR, >Z**I.VAR, and >Z**.COV (where ** is XX, XY, YX, or YY)
% keywords are provided for the real variance, imaginary variance, and
% covariance, respectively, for each tensor element. (Section 13.0).
%
% However, the usual practice is to calculate a simplified "variance" as
% defined by Gamble, 1978, pp. 66-72^*. This is a real number (an estimate of
% the average of the variances of the real and imaginary parts) that is the
% radius of a circle about the tensor element in the complex plane. This
% circle approximates the error ellipse and for most purposes is an
% adequate estimate of the statistical uncertainty. The >Z**.VAR keyword is
% provided (Section 13.0) for this parameter. This is typically the only
% error estimate provided for impedance tensor components.
%
% Gamble, 1978: https://escholarship.org/content/qt1qf5644g/qt1qf5644g_noSplash_db0b80f10299d43aac0d4bba46bb03bf.pdf#page=76