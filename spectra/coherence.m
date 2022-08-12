function [Cxy,fe] = coherence(x,y,opts,n)

if nargin < 4
    % To allow eval frequencies to be same as that for x having a different
    % number of time points.
    n = size(x,1);
end

N = size(x,1);
if mod(N,2) == 0
    Np = N/2 + 1;
else
    Np = (N-1)/2 + 1;    
end

ftx = fft(x);
ftx = ftx(1:Np,:);
fty = fft(y);
fty = fty(1:Np,:);

if nargin > 2
    winfn = opts.fd.window.function;
    [fe,Ic,Ne] = evalfreq(n,opts.fd.evalfreq.functionargs{:});
else
    winfn = @rectwin;
    Ic = 1:Np;
    Ne = zeros(Np,1);
    [~,fe] = fftfreq(N);
end    

Cxy(1,size(x,2)) = NaN;

for j = 2:length(Ic)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Caution - code below is duplicated in transferfnFD()
    W = winfn(2*Ne(j)+1);
    W = W/sum(W);
    r = [Ic(j)-Ne(j):Ic(j)+Ne(j)];
    % End duplicated code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for k = 1:size(x,2)
        sxx = sum(W.*(ftx(r,k).*conj(ftx(r,k))));
        syy = sum(W.*(fty(r,k).*conj(fty(r,k))));
        sxy = sum(abs(W.*(ftx(r,k).*conj(fty(r,k)))));
        Cxy(j,k) = sxy.^2./(sxx*syy);
    end

end