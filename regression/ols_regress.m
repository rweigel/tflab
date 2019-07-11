function [Z,W,stats] = ols_regress(ftB,ftE,varargin)

    Z = regress(ftE,ftB);
    W = [];
    stats = struct();

end
