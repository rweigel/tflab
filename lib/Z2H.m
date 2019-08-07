function [h,t] = z2h(Z,opts)
% Z2H Convert frequency domain transfer function to impulse response
%
%   [h,t] = Z2H(Z) returns h = ifft(Z) and t = N*fftfreq(N) where N is
%   number of elements in first non-singleton dimension.
% 
%   See also H2Z, ZINTERP.

if nargin < 2
    opts = struct();
end

assert(ndims(Z) <= 2,'Required: ndims(Z) <= 2');

if nargin == 1
    h = ifft(Z);
    if numel(Z) == length(Z) && size(Z,2) > 1
        N = size(Z,2);
        t = N*fftfreq(N);
    else
        N = size(Z,1);
        t = round(N*fftfreq(N))';
    end
    return;
end

% Compute impulse response
h = ifft(Zfull);
t = N*fftfreq(N)';

assert(isreal(h),...
            ['Computed impulse response has imaginary component.'...
            'Check calculation of Z']);

