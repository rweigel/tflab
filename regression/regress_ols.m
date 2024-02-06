function [Z,dZ,Info] = regress_ols(ftE,ftB,algorithm)
%REGRESS_OLS - Ordinary least-squares regression
%
%   Z = OLS_REGRESS(E, B) and Z = OLS_REGRESS(E, B, 'backslash')
%   returns the OLS solution to E = Z*B using Z = ftB\ftE. See HELP SLASH.
%
%   Z = OLS_REGRESS(E, B, 'regress') returns Z = regress(E, B); the
%   Statistics and Machine Learning Toolbox is required.
%
%   Z = OLS_REGRESS(E, B, 'analytic') returns Z = regress_ols_analytic(E, B)
%
%   Z = OLS_REGRESS(E, B, 'regress-real') returns Z computed using
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

W = ones(size(ftE,1),1);

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
    [Z,Info.ZCL95,Info.Residuals,Info.RINT,Info.STATS] = regress(ftE,ftB);
    if 1
        if offset
            Info.ZCL95 = Info.ZCL95(1:end-1,:);
        end
        % Z has is one column with rows of component. Last element is offset.
        % ZCL95 has rows of component. First column is lower limit, second is
        % upper limit. I don't think MATLAB is handling confidence limits
        % with complex Z properly. The imaginary part of the upper and lower
        % confidence limits on Z are the same as the imaginary part of Z.
        Info.ZCL95l = transpose(Info.ZCL95(:,1));
        Info.ZCL95u = transpose(Info.ZCL95(:,2));
        tmp = [Info.ZCL95l(1),Z(1),Info.ZCL95u(1)];
        fprintf('real(ZCLl(1), Z(1), ZCLu(1))')
        real(tmp)
        fprintf('imag(ZCLl(1), Z(1), ZCLu(1))')
        imag(tmp)
    else
        Info = rmfield(Info,'ZCL95');        
    end
end

if strcmp(algorithm,'regress-analytic')
    [Z,Info.Residuals] = regress_ols_analytic(ftE,ftB);
end

if offset
    dZ = Z(end);
    Z = Z(1:end-1);
else
    dZ = [];
end

Z = transpose(Z);

if strcmp(algorithm,'regress-real')
    % Calculation using only real numbers
    %   Ec = Zc*Bc for complex-valued scalars (subscript c) can be written as
    %
    %   Er + jEi = (Zr + jZi)*(Br + jBi),
    %
    %   where subscript r (i) is the real (imaginary) component. The two
    %   independent equations are
    %
    %   Er = Zr*Br - Zi*Bi
    %   Ei = Zr*Bi + Zi*Br,
    %
    %   or, in MATLAB matrix notation,
    %
    %   [Er ; Ei] = [Br , -Bi ; Bi , Br]*[Zr ; Zi]
    %  
    %   And so Zc = regress(Ec,Bc)
    %
    %   is equivalent to 
    % 
    %   Zc = Z(1:end/2,:) + sqrt(-1)*Z(end/2+1:end,:);
    %
    %   where
    %
    %   Z = regress([real(Ec) ; imag(Ec)],
    %               [real(Bc), -imag(Bc) ; imag(Bc), real(Bc)])
    [Zreal,~,Rreal] = regress(...
                        [real(ftE) ; imag(ftE)], ...
                        [real(ftB), -imag(ftB) ; imag(ftB), real(ftB)]);
    Z = Zreal(1:end/2,:) + sqrt(-1)*Zreal(end/2+1:end,:);
    R = Rreal(1:end/2,:) + sqrt(-1)*Rreal(end/2+1:end,:);
end
    