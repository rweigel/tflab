function varargout = trimnans(varargin)
%TRIMNANS Remove NaNs at start and end of matrices
%
%  [Ar,Br,...] =  TRIMNANS(A,B,...)
%
%  Example
%       A = [1 1 1; 2 2 2; 3 3 3; 4, 4, 4]
%       B = [nan; 99; 99; nan]
%       [a,b] = trimnans(A,B)
%
%       A =
%           1     1     1
%           2     2     2
%           3     3     3
%           4     4     4
%
%       B =
%           NaN
%           99
%           99
%           NaN
%
%       a =
%           2     2     2
%           3     3     3
%
%       b =
%           99
%           99
%
%  See also naninterp1.

nr = size(varargin{1},1);
for m = 2:length(varargin) % Loop over matrices
    assert(size(varargin{m},1) == nr,...
            'Number of rows of all matrices must be the same');
end

for m = 1:length(varargin) % Loop over matrices
    if length(size(varargin{m})) > 2
        error('Input matrix must have less than 3 dimensions');
    end
    for j = 1:size(varargin{m},2) % Loop over columns
        f = find(~isnan(varargin{m}(:,j)),1,'first');
        F{m}(j) = 1;
        if ~isempty(f)
            F{m}(j) = f;
        end
        l = find(~isnan(varargin{m}(:,j)),1,'last');
        L{m}(j) = nr;
        if ~isempty(l)
            L{m}(j) = l;
        end
    end
end

Fmax = max(cat(2,F{:}));
Lmin = min(cat(2,L{:}));

for m = 1:length(varargin)
    varargout{m} = varargin{m}(Fmax:Lmin,:);
end
