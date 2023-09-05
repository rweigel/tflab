function [X,W] = tdwindow(X,winfn,varargin)
%TDWINDOW Time domain window
%
%   [X,W] = TDWINDOW(X,winfn) applies winfn to each column of X:
%
%    W = winfn(size(X,1));
%    X = X.*repmat(W,1,size(X,2));
%
%   [X,W] = TDWINDOW(X,winfn,...) uses winfn arguments according to
%
%    W = winfn(size(X,1),varargin{:});
%    
%    See `help windows` for a list of window functions. See
%    `help dpss` for information on the dpss function.
%
%    See also TDWINDOW_DEMO, WINDOWS, DPSS.

if strcmp(winfn, 'rectwin')
    W = X;
    return;
end

if nargin < 3
    W = winfn(size(X,1));
else
    Ws = winfn(size(X,1),varargin{:});
    W = sum(Ws,2);
end

X = X.*repmat(W,1,size(X,2));
