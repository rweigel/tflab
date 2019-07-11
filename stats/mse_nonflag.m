function MSE = mse_nonflag(y,f,FLAG,COND)
%MSE_NONFLAG  Mean square error between non-FLAG elements
%
%   MSE(A,P) computes the mean square error for rows where both A and P
%   have non-NaN elements.
%
%   For matrices, MSE(A,P) returns a column vector with elements containing
%   the MSE between the respective columns of A and P.
%
%   MSE_NONFLAG(A,P,FLAG,COND) See IS_FLAG for a definition of FLAG and COND.
%
%   If all elements along dimension DIM are equal to FLAG, a fill of FLAG is
%   used.
%   
%   See also MSE, *_NONFLAG.

% R.S. Weigel, 04/02/2004.
  
if (nargin < 3)
  FLAG = NaN;
end
if (nargin < 4)
  COND = 1;
end

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

for i = 1:size(y,2)
  Iy = is_flag(y(:,i),FLAG,COND);
  If = is_flag(f(:,i),FLAG,COND);
  I  = find((Iy == 0) & (If == 0));
  if (length(I) > 1)
      MSE(1,i) = mse(y(I,i),f(I,i));
  else
      MSE(1,i) = FLAG;      
  end
end
