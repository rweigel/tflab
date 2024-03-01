function MSE = mse(y,f)
%MSE Mean square error between actual and predicted vector.
%
%   MSE(A,P) returns (1/N)*sum( (A-P).^2 )
%
%   For matrices, MSE(A,P) returns a column vector with elements containing
%   the MSE between the respective columns of A and P.
%
%   See also MSE_NONFLAG, ARV.

if (prod(size(y)) ~= prod(size(f)))
    error('Inputs must be the same size');
end

if size(y,1) == 1
  y = y';
  f = f';
  if (length(y) < 2)
    error('Inputs must have 2 or more elements');
  end
end

MSE = (sum((y - f).^2,1))/size(y,1);


