function [Z,W,stats] = ols_analytic(ftB,ftE,varargin)
    
    assert(size(ftB,2) <= 2,'ols_analytic() only one or two inputs');

    if size(ftB,2) == 1
        % Ex = ZxxBx
        Z = sum(ftE(:,1).*conj(ftB(:,1)))/sum(ftB(:,1).*conj(ftB(:,1)));
    else
        % Ex = ZxxBx + ZxyBy
        BxBx = sum(ftB(:,1).*conj(ftB(:,1))); 
        BxBy = sum(ftB(:,1).*conj(ftB(:,2)));

        ByBy = sum(ftB(:,2).*conj(ftB(:,2))); 
        ByBx = sum(ftB(:,2).*conj(ftB(:,1)));

        ExBx = sum(ftE(:,1).*conj(ftB(:,1))); 
        ExBy = sum(ftE(:,1).*conj(ftB(:,2))); 

        DET =  BxBy*ByBx - BxBx*ByBy;
        Zxx = (ExBy*ByBx - ExBx*ByBy)/DET;
        Zxy = (ExBx*BxBy - ExBy*BxBx)/DET;
        Z = [Zxx,Zxy];
    end
    W = [];
    stats = struct();    
end