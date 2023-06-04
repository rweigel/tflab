function f = fftfreqp(N, dt)
%FFTFREQP Positive Discrete Fourier Transform sample frequencies
%
%   f = FFTFREQP(N) returns
%
%   f = (1/N)*[0, 1, ..., N/2-1] if N is even
%   f = (1/N)*[0, 1, ..., (N-1)/2] if N is odd
%
%   f = FFTFREQP(N, dt) returns FFTFREQP(N)/dt, where dt is the time between
%   samples.
%   
%   See also FFTFREQ, EVALFREQ.

if nargin == 1
    dt = 1;
end

assert( N >= 0,'N must be >= 0');

if N == 0
    f = [];
    return;
end

if mod(N,2) == 0
    fp = 0:N/2-1;
else
    fp = 0:(N-1)/2;
end

f = fp'/(N*dt);
