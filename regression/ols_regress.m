function [Z,W,R] = ols_regress(ftB,ftE,varargin)
%OLS_REGRESS - Ordinary least-squares regression
%
%   Z = OLS_REGRESS(B, E) returns Z = regress(B, E);
%
%   Z = OLS_REGRESS(B, E, 1) returns Z = regress(B, E) computed using
%   
%   Zreal = regress([real(Ec) ; imag(Ec)], 
%                   [real(Bc), -imag(Bc) ; imag(Bc), real(Bc)])
%
%   and
%
%   Z = Zreal(1:end/2,:) + sqrt(-1)*Zreal(end/2+1:end,:);
%
%   Notes:
%
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

realvalued = 0;
loglevel = 0;
if nargin > 2
    opts = varargin{1};
    if isfield(opts,'realvalued')
        realvalued = opts.realvalued;
    end
    if isfield(opts,'loglevel')
        loglevel = opts.loglevel;
    end
end

sE = size(ftE);
sB = size(ftB);

if realvalued == 1
    % Should give same answer as realvalued == 0
    if loglevel > 0
        msg = 'Using regress(y, x) with real-valued matrices;';
        logmsg('%s size(y)=[%d,%d]; size(x) = [%d,%d]\n',...
                msg,sE(1),sE(2),sB(1),sB(2));
    end
    % Calculation using only real numbers
    [Zreal,~,Rreal] = regress([real(ftE) ; imag(ftE)], ...
                    [real(ftB), -imag(ftB) ; imag(ftB), real(ftB)]);
    Z = Zreal(1:end/2,:) + sqrt(-1)*Zreal(end/2+1:end,:);
    R = Rreal(1:end/2,:) + sqrt(-1)*Rreal(end/2+1:end,:);
else
    if loglevel > 0
        msg = 'Using regress(y, x) with complex-valued matrices;';
        logmsg('%s size(y)=[%d,%d]; size(x) = [%d,%d]\n',...
                msg,sE(1),sE(2),sB(1),sB(2));
    end
    [Z,~,R] = regress(ftE,ftB);
end

W = [];
    