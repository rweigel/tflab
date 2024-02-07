function [Z,dZ,Info] = regress_ols(ftE,ftB,algorithm)
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
%   See also REGRESS_TLS, REGRESS_ROBUSTFIT_TFLAB, and
%   REGRESS_ROBUSTFIT_MATLAB.

if nargin < 3
    algorithm = 'backslash';
end

algorithms = {'backslash','regress','regress-real','regress-analytic'};
emsg = sprintf('algorithm must be one of: %s',strjoin(algorithms,', '));
assert(any(strcmp(algorithm,algorithms)),emsg);

offset = 1;
if offset
    ftB = [ftB,ones(size(ftB,1),1)];
end

Info = struct();

if strcmp(algorithm,'backslash')
    Z = ftB\ftE;
    Info.Residuals = ftE - ftB*Z;
end

if strcmp(algorithm,'regress')
    % As demonstrated in regress_demo_complex.m, ZCL is wrong
    % when ftE or ftB are complex.
    [Z,ZCL95,Info.Residuals,Info.RINT,Info.STATS] = regress(ftE,ftB);
end

if strcmp(algorithm,'regress-analytic')
    [Z,Info.Residuals] = regress_ols_analytic(ftE,ftB);
end

if strcmp(algorithm,'regress-real')
    if isreal(ftE) && isreal(ftB)
        if offset
            ftB = ftB(:,1:end-1);
        end
        [Z,dZ,Info] = regress_ols(ftE,ftB,'regress');
        return;
    end
    % Calculation using only real numbers
    %   Ec = Zc*Bc for complex-valued scalars (subscript c) can be written as
    %
    %   Er + jEi = (Zr + jZi)*(Br + jBi),
    %
    %   where subscript r (i) indicates a real (imaginary) value and
    %   Er, Ei, Zr, Zi, Br, Bi are real numbers. Two independent equations are
    %
    %   Er = Zr*Br - Zi*Bi
    %   Ei = Zr*Bi + Zi*Br,
    %
    %   or, in MATLAB matrix notation, real matrices E and B are
    %
    %   E = [Er ; Ei]
    %   B = [Br , -Bi ; Bi , Br]*[Zr ; Zi]
    %  
    %   Zc = regress(Ec,Bc) and Z = regress(E,B) are related by
    % 
    %   Zc = Z(1:end/2,:) + sqrt(-1)*Z(end/2+1:end,:);\
    E = [real(ftE) ; imag(ftE)];
    B = [real(ftB), -imag(ftB) ; imag(ftB), real(ftB)];
    [Z,ZCL95,R,RINT,STATS] = regress(E,B);
    ZCL95 = ZCL95(1:end/2,:) + sqrt(-1)*ZCL95(end/2+1:end,:);
    Z = Z(1:end/2,:) + sqrt(-1)*Z(end/2+1:end,:);
    R = R(1:end/2,:) + sqrt(-1)*R(end/2+1:end,:);
end

if offset
    dZCL95 = ZCL95(end,:);
    ZCL95 = ZCL95(1:end-1,:);
    dZ = Z(end);
    Z = Z(1:end-1);
else
    dZ = [];
end

Info.ZCL95l = transpose(ZCL95(:,1));
Info.ZCL95u = transpose(ZCL95(:,2));
Info.dZCL95l = transpose(dZCL95(:,1));
Info.dZCL95u = transpose(dZCL95(:,2));

Z = transpose(Z);

    