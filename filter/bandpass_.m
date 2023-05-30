function [x,aib] = bandpass_(x,fb)
%BANDPASS_ - Bandpass filter signal in frequency domain
%
%  y = bandpass_(x,fb) Given a time series x and low/high bandpass
%  frequencies fb = [fl,fh], sets to zero the fourier coefficients of x 
%  outside of the range fb and returns the inverse fourier transform of
%  the modified coefficients.
%
%  See also BANDPASS, BANDPASS_TEST.

flip = 0;
if size(x,1) == 1
    flip = 1;
    x = x';
end

aib = fft(x);
N   = size(x,1);
f   = fftfreq(N);

for i = 1:size(x,2)
    if any(isnan(x(:,i)))
        warning('Column %d of x has one or more NaNs. Output for such columns be all NaNs',i);
    end
end

if length(fb) == 1
    Ib = find(fb == abs(f));
    if isempty(Ib)
        [~,Ib] = min(abs(fb - f));  
        warning('No exact match for requested frequency. Using nearest frequency = %f',f(Ib));
    end
else
  if fb(2) == fb(1)
    warning('fb(2) == fb(1). Calling bandpass with a single frequency.');
	[x,aib] = bandpass_(x,fb);
	return;
  end
  assert(fb(2) > fb(1),'fb(2) must be greater than fb(1)');
  Ib = find( abs(f) < fb(2) & abs(f) > fb(1) );
end

Io = ones(size(f));  % Io has index of omitted frequencies
Io(Ib) = 0;          % Mask band frequencies
aib(Io == 1,:) = 0;  % Set DFT coeficients to zero outside band

x = ifft(aib);
if flip
    x = x';
end
