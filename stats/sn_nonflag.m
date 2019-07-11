function sn = sn_nonflag(x,y)
%SN_NONFLAG
%
%  SN_NONFLAG(A,P) returns sum(x.^2)./sum( (x-y).^2 )
%

sn = sum(x.^2)./sum( (x-y).^2 );

