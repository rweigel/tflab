function rho = z2rho(f, Z)
%Z2RHO
%
%   rho = Z2RHO(f, Z) returns apparent resistivity in Ohm·m assuming
%   f in Hz and Z in (mV/km)/nT so that rho = |Z|^2/(5*f).
%
%   By definition,
%
%   rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
%
%   or, equivalently using E_i = Z_{ij}*B_j
%
%   rho_ij = (mu_o/omega)*|Z_ij|^2
%
%   With mu_o = 4*pi*e-7 N/A^2, omega = 2*pi*f, and
% 
%   (mV/km)/nT = 10^{3} (V/m)/T
%
%   rho_ij = |Z_ij|^2/(5*f) in Ohm·m
%
%   See also Egbert, Booker, and Schultz 1992 pg 15,116.

rho = abs(Z).^2./(5*f);
