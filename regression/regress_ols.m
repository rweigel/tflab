function [Z,Info] = regress_ols(ftE,ftB,algorithm)
%REGRESS_OLS - Ordinary least-squares regression
%
%   Z = REGRESS_OLS(E, B) and Z = REGRESS_OLS(E, B, 'backslash')
%   returns the OLS solution to E = Z*B using Z = ftB\ftE. See HELP SLASH.
%
%   Z = REGRESS_OLS(E, B, 'regress') returns Z = regress(E, B); the
%   Statistics and Machine Learning Toolbox is required.
%
%   Z = REGRESS_OLS(E, B, 'analytic') returns Z = regress_ols_analytic(E, B)
%
%   Z = REGRESS_OLS(E, B, 'regress-real') returns Z computed using
%   regress() and real values according to
%
%     Zreal = regress([real(Ec) ; imag(Ec)],
%                     [real(Bc), -imag(Bc) ; imag(Bc), real(Bc)])
%
%     and
%
%     Z = Zreal(1:end/2,:) + sqrt(-1)*Zreal(end/2+1:end,:);
%
%   See also REGRESS_TLS, REGRESS_ROBUSTFIT_TFLAB, REGRESS_ROBUSTFIT_MATLAB
%   REGRESS_DEMO, REGRESS_DEMO_COMPLEX.

if nargin < 3
    algorithm = 'backslash';
end

algorithms = {'backslash','regress','regress-real','regress-analytic'};
emsg = sprintf('algorithm must be one of: %s',strjoin(algorithms,', '));
assert(any(strcmp(algorithm,algorithms)),emsg);

Info = struct();

if strcmp(algorithm,'backslash')
    Z = ftB\ftE;
    Info.Residuals = ftE - ftB*Z;
end

if strcmp(algorithm,'regress')
    if ~isreal(ftE) || ~isreal(ftB)
        logmsg('ftIn or ftOut is complex. Using regress-real.\n');
        % As demonstrated in regress_demo_complex_bug.m, ZCL is wrong when
        % ftE or ftB is complex.
        [Z,Info] = regress_ols(ftE,ftB,'regress-real');
        return
    end
    [Z,ZCL95,Info.Residuals,Rint,Stats] = regress(ftE,ftB);
end

if strcmp(algorithm,'regress-analytic')
    [Z,Info.Residuals] = regress_ols_analytic(ftE,ftB);
end

if strcmp(algorithm,'regress-real')
    % Calculation using only real numbers
    %   Ec = Zc*Bc for complex-valued scalars (subscript c) can be written as
    %
    %   Er + jEi = (Zr + jZi)*(Br + jBi) + dZr + jdZi,
    %
    %   where subscript r (i) indicates a real (imaginary) value and
    %   Er, Ei, Zr, Zi, Br, Bi are real numbers. Two independent equations are
    %
    %   Er = Zr*Br - Zi*Bi + dZr
    %   Ei = Zr*Bi + Zi*Br + dZi
    %
    %   or, in matrix notation, real matrices E and B are
    %
    %   E = [Er]
    %       [Ei]
    %   B = [Br  -Bi  1  0][Zr]
    %       [Bi   Br  0  1][Zi]
    %                      [dZr]
    %                      [dZi]
    %   Zc = regress(Ec,Bc) and Z = regress(E,B) are related by
    %
    %   Zc = Z(1:end/2,:) + sqrt(-1)*Z(end/2+1:end,:);

    E = [real(ftE); ...
         imag(ftE)];
    B = [real(ftB), -imag(ftB);...
         imag(ftB),  real(ftB)];
    [Z,ZCL95,Residuals,Rint,Stats] = regress(E,B);
    Z = Z(1:end/2) + sqrt(-1)*Z(end/2+1:end);
    ZCL95 = ZCL95(1:end/2,:) + sqrt(-1)*ZCL95(end/2+1:end,:);
    Info.Residuals = Residuals(1:end/2) + sqrt(-1)*Residuals(end/2+1:end);
end

if exist('ZCL95','var')
    Info.ZCL95l = transpose(ZCL95(:,1));
    Info.ZCL95u = transpose(ZCL95(:,2));
end

Z = transpose(Z);
