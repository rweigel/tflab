function ARV = arv_nonflag(y,f,FLAG,COND)
%ARV_NONFLAG Average relative variance between non-FLAG elements
% 
%   ARV_NONFLAG(A,P) computes the ARV between non-NaN elements of A and P.
%
%   ARV_NONFLAG(A,P,FLAG,COND) See IS_FLAG for a definition of FLAG and COND.
%
%   See also ARV, MSE, *_NONFLAG.

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
    meany    = repmat(mean(y(I,i),1),length(I),1);
    N        = sum((y(I,i) - f(I,i)).^2,1);
    D        = sum((y(I,i) - meany).^2,1);
    ARV(1,i) = N./D;
  else
    ARV(1,i) = FLAG;      
  end
end
