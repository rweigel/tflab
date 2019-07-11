function [Z,W,stats] = robust_robustfit(ftB,ftE,varargin)

    [nr,nc] = size(ftE);
    if nr <= nc % Not enough data to use robustfit
        Z = regress(ftE,ftB);
        stats = [];
        W = [];
    else
        [Z,stats] = robustfit(ftB,ftE,varargin{:});
        W = stats.w;
    end
    
end