function I = is_flag(x,FLAG,COND)
%IS_FLAG Tests elements of array to determine if they qualify as flagged
%
%    I = IS_FLAG(x) Same as ISNAN(x)
%  
%    I = IS_FLAG(x,FLAG,COND)
%
%    Default FLAG is NaN.
%
%    COND = -1 condition is equivalent to ISNAN(x); FLAG must be [] or NaN
%    COND =  0 flag condition is x >= FLAG.
%    COND =  1  condition is x == FLAG (default).  
%    COND =  2  condition is abs(x) >= FLAG.
%
%    See also *_NONFLAG.

if (nargin < 2)
  FLAG = NaN;
  COND = 1;
end
if (nargin < 3)
  COND = 1;
end
if (isempty(FLAG))
   FLAG = NaN;
end

if (COND == -1) & (~isnan(FLAG))
   error('For COND = -1, FLAG must be [] or NaN');
end

I = zeros(size(x));

if (COND == -1)
  I = isnan(x);
  return;
end
if (COND == 0)
  I = (x >= FLAG);
  return;
end
if (COND == 1)
  if (isnan(FLAG))
    I = isnan(x);
  else
    I = (x == FLAG);
  end
  return;
end
if (COND == 2)
  I = (abs(x) >= FLAG);
  return;
end
