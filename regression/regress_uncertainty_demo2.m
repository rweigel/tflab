% Manually compupte error bars on regression coefficients for 1 and 2-d
% regression.

clear
Np = 100;
p = 1-.05/2;
sigma = 0.1;
bo = 1;

% 1-D
x = [1:Np]'/Np;
err = sigma*randn(Np,1);
y = bo + x + err;
X = [x,ones(Np,1)];
%err = err-mean(err);

if 1
    % 2-D
    x1 = randn(Np,1);
    x2 = randn(Np,1);
    err = sigma*randn(Np,1);
    y = bo + x1 + x2 + err;
    X = [x1,x2,ones(Np,1)];
end


[b,bint,r,rint,stats] = regress(y,X,0.05);
b
bint
stats

% Based on inspection of code from regress.m, 
% Draper and Smith, 1981, Applied Regression Analysis, p94 (which regress.m references)
% and https://www.mathworks.com/matlabcentral/fileexchange/202-mregress
% (but see 2009 comment by D'Errico in Reviews).

b = X\y

residuals = y - X*b;
[n, k] = size(X);
s2 = transpose(residuals) * residuals / (n - k) ;

XTXI = inv(X'*X);

tval = tinv(1-0.05/2, n-k);
se = sqrt(s2*diag(XTXI));

bintm = [b-tval*se, b+tval*se] % Manual calculation

bint - bintm


