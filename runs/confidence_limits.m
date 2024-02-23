function S = confidence_limits(S)

    if isfield(S,'Regression')
        if isfield(S.Regression.ErrorEstimates,'Parametric')
            tmp = S.Regression.Parametric.ErrorEstimates;
            S.Metrics.ErrorEstimates.Parametric = ...
                compute_derived_limits(S, tmp, 'Parametric');
        end
        if isfield(S.Regression.ErrorEstimates,'Bootstrap')
            tmp = S.Regression.Bootstrap.ErrorEstimates;
            S.Metrics.ErrorEstimates.Bootstrap = ...
                compute_derived_limits(S, tmp, 'Bootstrap');
        end
    end

end

function compute_derived_limits(S, ErrorEstimates)

if ~isfield(S,'ZVAR')
    %   |Z| = sqrt(Zr^2 + Zi^2)
    % standard propagation of errors is
    %   Δ|Z| = (|Zr|/|Z|)ΔZr + (|Zi|/|Z|)ΔZr
    %        = (1/|Z|)(|Zr|ΔZr + |Zi|ΔZr)
    % or
    %   Δ|Z| = sqrt( [(Zr/|Z|)ΔZr]^2 + [(|Zi|/|Z|)ΔZi]^2 )
    %   Δ|Z| = (1/|Z|)sqrt([ZrΔZr]^2 + [ZiΔZi]^2 )
    %
    % TODO: Derive ZCL68{u,l} from ZCL95 for parametric

    % Zmag = |Z|
    Zmag(:,comp) = abs(S.Z(:,comp));

    % Z95CL{l,u} are (complex-valued) error estimates on Z, not |Z|.
    ZCL95l(:,comp) = ErrorEstimates.ZCL95l(:,comp);
    ZCL95u(:,comp) = ErrorEstimates.ZCL95u(:,comp);

    delta_Zr = real(ZCL95u(:,comp) - ZCL95l(:,comp));
    delta_Zi = imag(ZCL95u(:,comp) - ZCL95l(:,comp));

    Zr = real(S.Z(:,comp));
    Zi = imag(S.Z(:,comp));

    %delta_Zmag(:,comp) = (1./Zmag).*(abs(Zr).*delta_Zr + abs(Zi).*delta_Zi);
    delta_Zmag(:,comp) = (1./Zmag).*sqrt((abs(Zr).*delta_Zr).^2 + (abs(Zi).*delta_Zi).^2);
    
    ZMAGCLl = Zmag(:,comp) - delta_Zmag/2;
    ZMAGCLu = Zmag(:,comp) + delta_Zmag/2;
else
    delta_Zmag = sqrt(ZVAR);
    ZVAR(:,comp) = S.ZVAR(:,comp);
    ZMAGCLl = S.Z(:,comp) - sqrt(2)*sqrt(ZVAR(:,comp));
    ZMAGCLu = S.Z(:,comp) + sqrt(2)*sqrt(ZVAR(:,comp));

    % ρ = |Z|^2/(5*f) => Δρ = 2|Z|ΔZ/(5*f)
    % See note below
    delta_Z = sqrt(2)*sqrt(ZVAR(:,comp));
    RHOVAR(:,comp) = 2*abs(S.Z(:,comp)).*delta_Z./(5*S.fe);
end

ErrorEstimates.ZMAGCLl = ZMAGCLl;
ErrorEstimates.ZMAGCLu = ZMAGCLu;

RHO(:,comp) = (abs(S.Z(:,comp)).^2)/(5*S.fe);
RHOLl(:,comp) = RHO(:,comp) - RHOVAR(:,comp)/2;
RHOLu(:,comp) = RHO(:,comp) - RHOVAR(:,comp)/2;

ErrorEstimates.RHO = RHO;
ErrorEstimates.RHOLl = RHOLl;
ErrorEstimates.RHOLu = RHOLu;

% Computing error bars from SPUD XML ZVAR.
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

end