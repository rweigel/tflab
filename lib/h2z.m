function [Z,f] = h2z(H)
%H2Z - Compute frequency domain transfer function from impulse response
%
%  [Z,f] = H2Z(H) returns Z = fft(H) at frequencies f given by
%  [~,f] = fftfreq(size(H,1)).
%
%  See also Z2H, FREQZ.

addpath([fileparts(mfilename('fullpath')),'/../fft']);

Z = fft(H);
N = size(H,1);

[~,f] = fftfreq(N);
Z = Z(1:length(f),:);
