function ARV = arv(y,f)
%ARV Average relative variance between actual and predicted vector.
% 
%   ARV(A,P) returns (1/N)*sum( (A-P).^2 )./(std(A,1))^2;
%
%   ARV is the mean square error scaled to the variance of the actual
%   vector, A (if A=P, ARV=0).  It can be interpreted as the fraction of the
%   variance of A that is not predicted by P.
%
%   For matrices, ARV(A,P) returns a column vector with elements containing
%   the ARV between the respective columns of A and P.
%   
%   See also ARV_NONFLAG, MSE.

% R.S. Weigel, 04/02/2004.
  
if (prod(size(y)) ~= prod(size(f)))
  error('Inputs must be the same size');
end
  
if (size(y,1) == 1)
  y = y';
  f = f';
  if (length(y) < 2)
    error('Inputs must have 2 or more elements');
  end
end

meany  = repmat(mean(y,1),size(y,1),1);
N      = sum((y - f).^2,1);
D      = sum((y - meany).^2,1);
ARV    = N./D;
