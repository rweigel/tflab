function [f,fu] = fftfreq(N, s)
%FFTFREQ Discrete Fourier Transform sample frequencies
%
%   f = FFTFREQ(N) returns
%
%   f = (1/N)*[0, 1, ..., N/2-1, -N/2, ..., -1] if N is even
%   f = (1/N)*[0, 1, ..., (N-1)/2, -(N-1)/2, ..., -1] if N is odd
%
%   f = FFTFREQ(N, s) returns FFTFREQP(N)/d, where d is interpreted
%   as the sampling frequency.
%
%   [f,fu] = FFTFREQ(...) returns the unique DFT frequencies:
%
%   fu = (1/N)*[0,1,...,N/2] if N is even
%   fu = (1/N)*[0,1,...,(N-1)/2] if N is even
%
%   (For N even, the highest frequency magnitude is 0.5 and negative by
%   convention. fu uses +0.5.)
%   
%   See also FFTFREQP.

if nargin == 1
    s = 1;
end

assert(N >= 0,'N must be >= 0');

if N == 0
    f = [];
    return;
end

if mod(N,2) == 0 % Even
    fp = 0:N/2-1;
    fn = -N/2:-1;
    fu = [fp,N/2]/(N*s);
else
    fp = 0:(N-1)/2;
    fn = -(N-1)/2:-1;
    fu = fp/(N*s);
end

f = [fp,fn]/(N*s);    
