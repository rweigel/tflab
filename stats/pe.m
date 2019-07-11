function PE = pe(y,f)
%PE Prediction efficiency
% 
%   PE(A,P) returns 1-ARV(A,P).
%
%   A PE can be interpreted as the fraction of the variance in A that is
%   predicted by P (if A == P, PE = 1). It can be interpreted as the
%   fraction of the variance of A that is predicted by P.
%
%   For matrices, PE(A,P) returns a column vector with elements containing
%   the ARV between the respective columns of A and P.
%   
%   See also ARV, MSE, CC.

PE = 1-arv(y,f); 
