function [Z,W,stats] = ols_analytic(ftB,ftE,varargin)
%OLS_ANALYTIC - TF estimate using analytic ordinary least-squares regression
% 
%   Z = OLS_ANALYTIC(B,E) uses an analytic formula to compute the OLS
%   solution for Z in the linear model
% 
%   E(:,1) = Z*B(:,1)   or   E(:,1) = Z*B(:,1:2)
%
%   where E and B are real or complex. Solutions are numerically equivalent
%   to OLS_REGRESS(B,E).
%
%   See also OLS_REGRESS.

assert(size(ftE,1) == 1,'Number of input columns must be 1.');
assert(size(ftB,2) <= 2,'Number of output columns must be 1 or 2');
assert(size(ftE,1) == size(ftB,1)','E and B must have the same number of rows');

if size(ftB,2) == 1
    % Ex = ZxxBx
    Z = sum(ftE(:,1).*conj(ftB(:,1)))/sum(ftB(:,1).*conj(ftB(:,1)));
else
    % Ex = ZxxBx + ZxyBy
    BxBx = sum(ftB(:,1).*conj(ftB(:,1))); 
    BxBy = sum(ftB(:,1).*conj(ftB(:,2)));

    ByBy = sum(ftB(:,2).*conj(ftB(:,2))); 
    ByBx = sum(ftB(:,2).*conj(ftB(:,1)));

    ExBx = sum(ftE(:,1).*conj(ftB(:,1))); 
    ExBy = sum(ftE(:,1).*conj(ftB(:,2))); 

    DET =  BxBy*ByBx - BxBx*ByBy;
    Zxx = (ExBy*ByBx - ExBx*ByBy)/DET;
    Zxy = (ExBx*BxBy - ExBy*BxBx)/DET;
    Z = [Zxx,Zxy];
end
W = [];
stats = struct();    
