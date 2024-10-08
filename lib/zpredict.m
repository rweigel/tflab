function [Ep,Z] = zpredict(Z,B,Nout)
%ZPREDICT Predict output given input and freq. domain transfer function
%
%  E = ZPREDICT(Z,B) uses the frequency domain transfer function Z
%  and the time series B and returns the output time series E.
%
%  If the number of columns of Z is larger than that of B, E will have
%  size(Z,2)/size(B,2) columns. The first set of size(B,2) columns in Z are
%  used with B to form the first column of E. The second set of size(B,2)
%  columns in Z are convolved with B to form the second column of E, etc.
%
%  See also ZPREDICT_TEST.

assert(size(B,1) == size(Z,1),'Required: size(B,1) == size(Z,1)');

if size(Z,2)/size(B,2) == Nout
    const_term = 0;
    Nin  = size(B,2);    
else
    if size(Z,2)/(1+size(B,2)) ~= Nout
        error('size(Z,2)/size(B,2) or size(Z,2)/(1+size(B,2)) must be equal to Nout = %d',Nout);
    end
    const_term = 1;
    Nin  = 1 + size(B,2);    
end

for j = 1:Nout
    zcols = (1:Nin) + (j-1)*Nin;
    if const_term == 1
        Ep(:,j) = sum(ifft( fft(B).*Z(:,zcols(1:end-1)) + Z(:,zcols(end)) ),2);
        continue;
    else
        Ep(:,j) = sum(ifft( fft(B).*Z(:,zcols) ),2);
    end
end

c = max(max(abs(imag(Ep))));
if max(abs(imag(Ep))) > eps
    msg = 'Computed prediction has imaginary component with max(abs(imag(H)))/eps = %.1f. Check calculation of Z.';
    warning(msg,c/eps);
end
Ep = real(Ep);
