function [dftu,fu] = fftu(x)

[~,fu] = fftfreq(size(x,1)); % Unique DFT frequencies

% # of unique frequency values.
Nu = length(fu);

dft = fft(x);
dftu = dft(1:Nu,:);
