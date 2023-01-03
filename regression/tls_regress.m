function [Z,W,R] = tls_regress(ftB,ftE,varargin)

if size(ftB,1) <= size(ftB,2)
    warning('Not enough data to use pca(). Using regress().');
    [Z,~,R] = regress(ftE,ftB);
else
    coeff = pca([ftB,ftE]);
    Z = -coeff(1:end-1,end)/coeff(end,end);
    R = ftE - ftB*Z;
end

W = [];