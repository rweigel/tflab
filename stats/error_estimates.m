function Metrics = error_estimates(S)

%logmsg('Computing additional error estimates.\n');
if isfield_(S,'Regression.ErrorEstimates.Parametric')
    ParametricErrorEstimates = S.Regression.ErrorEstimates.Parametric;    
end
if isfield_(S,'Regression.ErrorEstimates.Bootstrap')
    Metrics.ErrorEstimates.Bootstrap = S.Regression.ErrorEstimates.Bootstrap;
end

for comp = 1:size(S.Z,2)

    if isfield(S,'ZVAR') && ~isfield(ParametricErrorEstimates,'ZCL95l')
        % This will occcur when S is based on content of EDI file, so no
        % TFLab regression information is available.
    
        logmsg('Found ZVAR field at top-level and no ZCL95l in Regression.ErrorEstimates.\n');
        logmsg('Computing ZMAGCL95{u,l} using it.\n');
        
        % The following assumes ZVAR is variance on |Z| estimate.
        ZMAGCL95l = abs(S.Z(:,comp)) - sqrt(S.ZVAR(:,comp));
        ZMAGCL95u = abs(S.Z(:,comp)) + sqrt(S.ZVAR(:,comp));

        sf = (1+1j)/sqrt(2);
        ZCL95l(:,comp) = S.Z(:,comp) - sf*sqrt(S.ZVAR(:,comp));
        ZCL95u(:,comp) = S.Z(:,comp) + sf*sqrt(S.ZVAR(:,comp));
        ZVAR(:,comp) = S.ZVAR(:,comp);
        
    else
        logmsg('Computing ZMAGCL{l,u} using ZCL95{u,l}.\n');
    
        %   |Z| = sqrt(Zr^2 + Zi^2)
        % standard propagation of errors is either
        %   Δ|Z| = (|Zr|/|Z|)ΔZr + (|Zi|/|Z|)ΔZr
        %        = (1/|Z|)(|Zr|ΔZr + |Zi|ΔZr)
        %
        % or the more conservative
        %
        %   Δ|Z| = sqrt( [(Zr/|Z|)ΔZr]^2 + [(|Zi|/|Z|)ΔZi]^2 )
        %   Δ|Z| = (1/|Z|)sqrt([ZrΔZr]^2 + [ZiΔZi]^2 )
        %
        % TODO: Derive ZCL68{u,l} from ZCL95 for parametric
            
        % Z95CL{l,u} are (complex-valued) error estimates on Z, not |Z|.
        ZCL95l(:,comp) = ParametricErrorEstimates.ZCL95l(:,comp);
        ZCL95u(:,comp) = ParametricErrorEstimates.ZCL95u(:,comp);
    
        delta_Zr_95 = real(ZCL95u(:,comp) - ZCL95l(:,comp));
        delta_Zi_95 = imag(ZCL95u(:,comp) - ZCL95l(:,comp));
    
        % Zmag = |Z|
        Zmag = abs(S.Z(:,comp));

        Zr = real(S.Z(:,comp));
        Zi = imag(S.Z(:,comp));
    
        sumsq = (abs(Zr).*delta_Zr_95).^2 + (abs(Zi).*delta_Zi_95).^2;
        delta_Zmag_95 = (1./Zmag).*sqrt(sumsq);
        
        ZMAGCL95l(:,comp) = Zmag - delta_Zmag_95/2;
        ZMAGCL95u(:,comp) = Zmag + delta_Zmag_95/2;
    
        % delta_ZMAG_95 is 2σ, so σ = delta_ZMAG_95/2.
        ZVAR(:,comp) = (delta_Zmag_95/2).^2;

        % φ = (180/pi) * atan(Zi/Zr) => 
        % Δφ = (180/pi) * (1/|Z|^2) sqrt( (ZiΔZr)^2 + (ZrΔZri)^2 )
        
        PHI = (180/pi) * atan2(imag(S.Z(:,comp)),real(S.Z(:,comp)));
        delta_PHI_95 = (180/pi)*(1./Zmag.^2)...
                      .*sqrt( (Zi.*delta_Zr_95).^2 + (Zr.*delta_Zi_95).^2);
        PHICL95l(:,comp) = PHI - delta_PHI_95/2;
        PHICL95u(:,comp) = PHI + delta_PHI_95/2;
    end

    % ρ = |Z|^2/(5*f) => Δρ = 2|Z|ΔZ/(5*f)
    delta_Z_95 = 2*sqrt(ZVAR(:,comp)); % See note 1. below.
    delta_RHO_95 = 2*abs(S.Z(:,comp)).*delta_Z_95./(5*S.fe);
    
    RHO = (abs(S.Z(:,comp)).^2)./(5*S.fe);
    RHOCL95l(:,comp) = RHO - delta_RHO_95/2;
    RHOCL95u(:,comp) = RHO + delta_RHO_95/2;


end

Metrics.ErrorEstimates.Parametric.ZVAR = ZVAR;
Metrics.ErrorEstimates.Parametric.ZCL95l = ZCL95l;
Metrics.ErrorEstimates.Parametric.ZCL95u = ZCL95u;
Metrics.ErrorEstimates.Parametric.ZMAGCL95l = ZMAGCL95l;
Metrics.ErrorEstimates.Parametric.ZMAGCL95u = ZMAGCL95u;
Metrics.ErrorEstimates.Parametric.RHOCL95l = RHOCL95l;
Metrics.ErrorEstimates.Parametric.RHOCL95u = RHOCL95u;
Metrics.ErrorEstimates.Parametric.PHICL95l = PHICL95l;
Metrics.ErrorEstimates.Parametric.PHICL95u = PHICL95u;

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
