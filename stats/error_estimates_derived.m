function Parametric = error_estimates_derived(fe,Z,ZCL95l,ZCL95u)

for comp = 1:size(Z,2)
    delta_Zr_95 = real(ZCL95u(:,comp) - ZCL95l(:,comp));
    delta_Zi_95 = imag(ZCL95u(:,comp) - ZCL95l(:,comp));

    % Zmag = |Z|
    Zmag = abs(Z(:,comp));

    Zr = real(Z(:,comp));
    Zi = imag(Z(:,comp));

    % If
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
    sumsq = (abs(Zr).*delta_Zr_95).^2 + (abs(Zi).*delta_Zi_95).^2;
    delta_Zmag_95 = (1./Zmag).*sqrt(sumsq);

    ZMAGCL95l(:,comp) = Zmag - delta_Zmag_95/2;
    ZMAGCL95u(:,comp) = Zmag + delta_Zmag_95/2;

    % delta_ZMAG_95 is 2σ, so σ = delta_ZMAG_95/2.
    % TODO: See note 2.; this is (probably) technically not ZVAR as
    %       used in EDI files.
    ZVAR(:,comp) = (delta_Zmag_95/2).^2;

    % φ = (180/pi) * atan(Zi/Zr) =>
    % Δφ = (180/pi) * (1/|Z|^2) sqrt( (ZiΔZr)^2 + (ZrΔZri)^2 )
    PHI = (180/pi) * atan2(imag(Z(:,comp)),real(Z(:,comp)));
    tmp = sqrt( (Zi.*delta_Zr_95).^2 + (Zr.*delta_Zi_95).^2);
    delta_PHI_95 = (180/pi)*(1./Zmag.^2).*tmp;
    PHICL95l(:,comp) = PHI - delta_PHI_95/2;
    PHICL95u(:,comp) = PHI + delta_PHI_95/2;

    % ρ = |Z|^2/(5*f) => Δρ = 2|Z|ΔZ/(5*f)
    delta_RHO_95 = 2*abs(Z(:,comp)).*delta_Zmag_95./(5*fe);
    RHO = (abs(Z(:,comp)).^2)./(5*fe);
    RHOCL95l(:,comp) = RHO - delta_RHO_95/2;
    RHOCL95u(:,comp) = RHO + delta_RHO_95/2;
end

Parametric.ZVAR = ZVAR;
Parametric.ZCL95l = ZCL95l;
Parametric.ZCL95u = ZCL95u;
Parametric.ZMAGCL95l = ZMAGCL95l;
Parametric.ZMAGCL95u = ZMAGCL95u;
Parametric.RHOCL95l = RHOCL95l;
Parametric.RHOCL95u = RHOCL95u;
Parametric.PHICL95l = PHICL95l;
Parametric.PHICL95u = PHICL95u;
