function [Z,W,stats] = robust_robustfit(ftB,ftE,varargin)

[nr,nc] = size(ftE);

if nr <= nc
    warning('Not enough data to use robustfit(). Using regress().');
    % TODO: Document why constraint is not nr < nc. The constraint nr <= nc
    % is from statrobustfit. One constraint is that nr > 1 so that a
    % standard deviation can be computed. But need to document why in
    % general nr <= nc is the appropriate constraint.
    Z = regress(ftE,ftB);
    stats = [];
    W = [];
else
    [Z,stats] = robustfit(ftB,ftE,varargin{:});
    W = stats.w;
end
