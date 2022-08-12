function [S,D,f] = psd(x,opts,n)

if nargin < 3
    % If time series x was trimmed, we can pass untrimmed length of x
    % so that the frequency grid is the same as that for the untrimmed x.
    n = size(x,1); 
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
    return
end

if nargin > 2
    winfn = opts.fd.window.function;
else
    winfn = @parzenwin;
end

[f,Ic,Ne] = evalfreq(n,opts.fd.evalfreq.functionargs{:});

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
