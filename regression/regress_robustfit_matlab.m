function [Z,R,W,stats] = regress_robustfit_matlab(ftE,ftB,varargin)

[nr,nc] = size(ftE);

if nr <= nc
    warning('Not enough data to use robustfit(). Using regress().');
    % TODO: Document why constraint is not nr < nc. The constraint nr <= nc
    % is from statrobustfit. One constraint is that nr > 1 so that a
    % standard deviation can be computed. But need to document why in
    % general nr <= nc is the appropriate constraint.
    [Z,~,R] = regress(ftE,ftB);
    stats = [];
    W = [];
else
    % Note robustfit() uses robustfit(input, output)
    % and regress() uses regress(output, input)
    [Z,stats] = robustfit(ftB,ftE,varargin{:});
    W = stats.w;
    % Unweighted residuals
    R = ftE - ftB*Z;
end

