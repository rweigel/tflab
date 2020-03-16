function [S,fe] = smoothSpectra(x,opts,n)

if nargin < 3
    % To allow eval frequencies to be same as that for x having a different
    % number of time points.
    n = size(x,1); 
end

if nargin > 2
    winfn = opts.fd.window.function;
else
    winfn = @parzenwin;
end

[fe,Ic,Ne] = evalfreq(n,opts.fd.evalfreq.functionargs{:});
 
N = size(x,1);
if mod(N,2) == 0
    Np = N/2 + 1;
else
    Np = (N-1)/2 + 1;    
end

ftB = fft(x);
ftB = ftB(1:Np,:);

for j = 2:length(Ic)

    W = winfn(2*Ne(j)+1);
    W = W/sum(W);
    r = Ic(j)-Ne(j):Ic(j)+Ne(j);

    for k = 1:size(ftB,2)
        S(j,k) = sum(W.*(ftB(r,k).*conj(ftB(r,k))));
    end

end
