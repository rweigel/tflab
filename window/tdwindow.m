function [X,W] = tdwindow(X,winfn,varargin)
%TDWINDOW Time domain window
%
%   [X,W] = TDWINDOW(X,winfn) computes
%
%    W = winfn(size(X,1));
%    X = X.*repmat(W,1,size(X,2));
%
%    See `help windows` for a list of window functions. See
%    `help dpss` for information on the dpss function.
%
%    See also windows, dpss.

% https://ccrma.stanford.edu/~jos/sasp/Kaiser_DPSS_Windows_Compared.html

if nargin < 3
    W = winfn(size(X,1));
else
    Ws = winfn(size(X,1),varargin{:});
    W = sum(Ws,2);
end

X = X.*repmat(W,1,size(X,2));
