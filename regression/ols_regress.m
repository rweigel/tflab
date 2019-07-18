function [Z,W,stats] = ols_regress(ftB,ftE,varargin)
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
%   Ec = Zc*Bc for complex-valued matrices (subscript c) can be written as
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

if realvalued == 1
    if loglevel > 0
        fprintf(['ols_regress(): Using regress() '...
                'with real-valued matrices\n']);
    end
    % Calculation using only reals numbers
    Zreal = regress([real(ftE) ; imag(ftE)], ...
                    [real(ftB), -imag(ftB) ; imag(ftB), real(ftB)]);
    Z = Zreal(1:end/2,:) + sqrt(-1)*Zreal(end/2+1:end,:);
else
    if loglevel > 0
        fprintf(['ols_regress(): Using regress() '...
                 'with complex-valued matrices\n']);
    end
    Z = regress(ftE,ftB);
end

W = [];
stats = struct();
    