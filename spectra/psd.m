function [S,D,f] = psd(x,opts)
%PSD Power spectra, Fourier coefficients, and frequencies
%
%  [S,D,f] = psd(x,opts).
%
%  Computes the power spectra (not density) S and Fourier coefficients
%  of signal x with N rows for positive frequencies according to
%
%    S = fft(x).*conj(fft(x))/(N/2)^2
%  
%  At a given frequency, a sinusoid with amplitude A and arbitrary phase
%  will have S = A^2 at that frequency.
%  

flip = 0;
if size(x,1) == 1
    flip = 1;
    x = transpose(x);
end

N = size(x,1);
if mod(N,2) == 0
    Np = N/2 + 1; % Keeps f = 0.5.
else
    Np = (N-1)/2 + 1;    
end

ftx = fft(x);
ftx = ftx(1:Np,:);

if nargin == 1 % Compute raw psd and fft
    S = ftx.*conj(ftx)/(N/2)^2;
    D = ftx;
    [~,f] = fftfreq(N);
    [S,D,f] = unflip(S,D,f,flip);
    return;
end

if nargin > 2
    winfn = opts.fd.window.function;
else
    winfn = @parzenwin;
end

[f,Ic,Ne] = evalfreq(N,opts.fd.evalfreq.functionargs{:});

%S(1,size(x,2)) = NaN;

for j = 1:length(Ic)

    W = winfn(2*Ne(j)+1);
    W = W/sum(W);
    r = Ic(j)-Ne(j):Ic(j)+Ne(j);
    for k = 1:size(ftx,2)
        S(j,k) = sum(W.*(ftx(r,k).*conj(ftx(r,k))));
        D(j,k) = sum(W.*(ftx(r,k)));
    end

end

S = S/(N/2)^2;
[S,D,f] = unflip(S,D,f,flip);

function [S,D,f] = unflip(S,D,f,flip)
    if flip == 1
        S = S';
        D = transpose(D);
        f = f';
    end
end
end
