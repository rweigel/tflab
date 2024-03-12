function Bootstrap = error_estimates_bootstrap(fe,Zb,Zm)
% Bootstrap confidence limits for Z
% Zb is a matrix of size Nb x size(Z,2)
% Returns a structure with fields ZCL95l, ZCL95u, and ZVAR
% ZCL95l and ZCL95u are estimates of the 95% confidence limits for Z
% ZVAR is an estimate of the variance of Z

% Consider also pivot CI:
% https://www2.stat.duke.edu/~banks/111-lectures.dir/lect13.pdf
for c = 1:size(Zb,2) % Columns are the transfer function components
    Bootstrap.ZCL95l(c) = prctile(real(Zb(:,c)),2.5) + 1j*prctile(imag(Zb(:,c)),2.5);
    Bootstrap.ZCL95u(c) = prctile(real(Zb(:,c)),97.5) + 1j*prctile(imag(Zb(:,c)),97.5);
    Bootstrap.ZVAR(c) = var(abs(Zb(:,c)),0,1);
    Bootstrap.ZMAGCL95l(c) = prctile(abs(Zb(:,c)),2.5);
    Bootstrap.ZMAGCL95u(c) = prctile(abs(Zb(:,c)),97.5);
    RHO = (abs(Zb(:,c)).^2)/(5*fe);
    Bootstrap.RHOCL95l(:,c) = prctile(RHO,2.5);
    Bootstrap.RHOCL95u(:,c) = prctile(RHO,97.5);
end

% Can't easily compute boostrap confidence limits for PHI b/c of wrapping.
% So compute it using error propagation using bootstrap ZCL95l and ZCL95u.
Parametric = error_estimates_derived(fe,Zm,Bootstrap.ZCL95l,Bootstrap.ZCL95u);
Bootstrap.PHICL95l = Parametric.PHICL95l;
Bootstrap.PHICL95u = Parametric.PHICL95u;
