function [X,W] = tdwindow(X,winfn)
%TDWINDOW
%
%   [X,W] = tdwindow(X,winfn) computes
%
%      W = winfn(size(X,1));
%      X = X.*repmat(W,1,size(X,2));

W = winfn(size(X,1));
X = X.*repmat(W,1,size(X,2));

    