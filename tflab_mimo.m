function S = tflab_mimo(B,E,t,opts)
% TFLAB_MIMO Frequency domain MIMO transfer function estimate
%
%  If B has Nb columns and E has Ne columns, Z will have 2*Nb*Ne columns,
%  with columns 1:Nb corresponding to the transfer function that gives
%  E(:,1) given B, columns Nb+1:2*Nb corresponding to the transfer function
%  that gives E(:,2) given B, etc.
%
%  In more familar notation, if the input and outputs are
%
%   Bx(t) = B(:,1)      By(t) = B(:,2)
%   Ex(t) = E(:,1)      Ey(t) = E(:,2)
%
%  the model equations are
%
%   Ex(f) = Zxx(f)Bx(f) + Zxy(f)By(f)
%   Ey(f) = Zyx(f)Bx(f) + Zyy(f)By(f)
%
%  and S = tflab_mimo(B,E) returns a structure S with a field Z such that
%
%   Zxx = Z(:,1)        Zxy = Z(:,2)
%   Zyx = Z(:,3)        Zyy = Z(:,4)
% 
%  with the rows of Z being estimates of the transfer function at the
%  evaluation frequencies given by field fe in S.
%
%  opts.td.window.width and opts.td.window.shift are ignored.

opts.td.window.width = NaN;
opts.td.window.shift = NaN;

S = tflab(B,E,t,opts);

