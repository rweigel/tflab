function [Ep,Z] = zpredict(Z,B)
%ZPREDICT Predict output given input and freq. domain transfer function
%
%  [E,Z] = ZPREDICT(Z,B) uses the frequency domain transfer function Z
%  and the time series B and returns the output time series E.
%  
%  If the number of columns of Z is larger than that of B, E will have
%  size(Z,2)/size(B,2) columns. The first set of size(B,2) columns in Z are
%  used with B to form the first column of E. The second set of size(B,2)
%  columns in Z are convolved with B to form the second column of E, etc.
%
%  See also ZPREDICT_DEMO.

assert(size(B,1) == size(Z,1),'Requires: size(B,1) == size(Z,1)');

Nout = size(Z,2)/size(B,2); % Number of output columns in E
Nin  = size(B,2);

assert(mod(size(Z,2),Nin) == 0,'size(Z,2)/size(B,2) must be an integer');

for j = 1:Nout
    for i = 1:size(B,2)
        c = (j-1)*size(B,2) + i;
        if i > 1
            Ep(:,j) = Ep(:,j) + ifft(fft(B(:,i)).*Z(:,c));
        else
            Ep(:,j) = ifft(fft(B(:,i)).*Z(:,c));            
        end
    end
end

if max(abs(imag(Ep))) > eps
    assert(isreal(Ep),...
            ['Computed impulse response has imaginary component.'...
            'Check calculation of Z']);
end
Ep = real(Ep);
