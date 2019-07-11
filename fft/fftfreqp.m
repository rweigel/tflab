function f = fftfreqp(N, s)
%FFTFREQ Positive Discrete Fourier Transform sample frequencies
%
%   f = FFTFREQP(N) returns
%
%   f = (1/N)*[0, 1, ..., N/2-1] if N is even
%   f = (1/N)*[0, 1, ..., (N-1)/2] if N is odd
%
%   f = FFTFREQP(N, d) returns FFTFREQP(N)/d, where d is interpreted
%   as the sampling frequency.
%   
%   See also FFTFREQ.

if nargin == 1
    s = 1;
end

assert(N>=0,'N must be >= 0');

if N == 0
    f = [];
    return;
end

if mod(N,2) == 0
    fp = 0:N/2-1;
else
    fp = 0:(N-1)/2;
end

f = fp/(N*s);
