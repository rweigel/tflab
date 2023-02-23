function [h,t] = z2h(Z)
% Z2H Convert frequency domain transfer function to impulse response
%
%   [h,t] = Z2H(z) returns h = ifft(z) and t = N*fftfreq(N) where N is
%   number of elements in first non-singleton dimension.
% 
%   See also H2Z, ZINTERP.

assert(ndims(Z) <= 2,'Required: ndims(Z) <= 2');

Inan1 = find(isnan(Z(1,:)));
if ~isempty(Inan1)
    warning('Z has NaN(s) at zero frequency. Setting to zero');
    Z(1,Inan1) = 0;
end
assert(all(~isnan(Z(:))),'Z has NaN values at non-zero frequencies');

h = ifft(Z);
if numel(Z) == length(Z) && size(Z,2) > 1
    N = size(Z,2);
    t = N*fftfreq(N);
else
    N = size(Z,1);
    t = round(N*fftfreq(N))';
end

assert(isreal(h),...
            ['Computed impulse response has imaginary component.'...
            'Check calculation of Z']);
