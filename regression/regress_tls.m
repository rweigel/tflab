function [Z,R] = regress_tls(Y,X)
%REGRESS_TLS - Total (or Orthogonal) Least Squares
%
%   [Z,R] = REGRESS_TLS(Y,X) compute the total least squares regression
%   of Y = Z*X using
%
%   coeff = pca([Y,X]);
%   Z = -coeff(1:end-1,end)/coeff(end,end);
%   R = Y - Z*X;
%
%   See also REGRESS_OLS.

if size(X,1) <= size(X,2)
    warning('size(Y,1) <= size(Y,2); not enough data to use pca(). Using regress().');
    [Z,~,R] = regress(Y,X);
else
    coeff = pca([X,Y]);
    Z = -coeff(1:end-1,end)/coeff(end,end);
    R = Y - X*Z;
end
