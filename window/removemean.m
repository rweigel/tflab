function varargout = removemean(varargin)
%REMOVEMEAN Remove column mean of non-nan elements
%
%  [Ar,Br,...] =  REMOVEMEAN(A,B,...)
%
% 
  
for c = 1:length(varargin)
    if length(size(varargin{c})) > 3
        error('Input matrix must have less than 4 dimensions');
    end
    varargout{c} = varargin{c}; % Preallocate
    for k = 1:size(varargin{c},3)
        for j = 1:size(varargin{c},2)
            Ig = find(~isnan(varargin{c}(:,j,k)));
            varargout{c}(:,j,k) = varargin{c}(:,j,k)-mean(varargin{c}(Ig,j,k));
        end
    end
end
