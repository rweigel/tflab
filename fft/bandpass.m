function [x,aib] = bandpass(x,fb)
%BANDPASS
%
%  y = bandpass(x,f)
%

aib = fft(x);
N   = length(x);
f   = fftfreq(N);

if (length(fb) == 1)
    Ib = find(fb == abs(f));
    if isempty(Ib)
        [~,Ib] = min(abs(fb - f));  
        warning('No exact match for requested period. Using nearest frequency = %f',f(Ib));
    end
else
  if fb(2) == fb(1)
    warning('fb(2) == fb(1). Calling bandpass with a single frequency.');
	[x,aib] = bandpass(x,fb);
	return;
  end
  assert(fb(2) > fb(1),'fb(2) must be greater than fb(1)');
  Ib = find( abs(f) < fb(2) & abs(f) > fb(1) );
end

Io = ones(size(f)); % Io has index of omitted frequencies
Io(Ib) = 0;         % Mask band frequencies
aib(Io == 1) = 0;   % Set DFT coeficients to zero outside band

x = ifft(aib);
