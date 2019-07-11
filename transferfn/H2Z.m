function [Z,f] = H2Z(H)
%H2Z
%
%  [Z,f] = H2Z(h) returns Z for frequencies in range fu, where
%  [~,fu] = fftfreq(size(h,1)).
%
%  See also Z2H.

Z = fft(H);
N = size(H,1);

[~,f] = fftfreq(N);
Z = Z(1:length(f),:);
