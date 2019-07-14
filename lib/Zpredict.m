function [Ep,Z] = Zpredict(Z,B,fe,opts)
%ZPREDICT Predict output given input and freq. domain transfer function
%
%  [Ep,Z] = Zpredict(Z,B) uses the frequency domain transfer function Z
%  and the time series B and returns the output timeseries Ep.
%  
%  If the number of columns of Z is larger than that of B, Ep will have
%  size(Z,2)/size(B,2) columns. The first set of size(B,2) columns in Z
%  are convolved with B to form the first column of Ep. The second set of
%  size(B,2) columns in Z are convolved with B to form the second column of
%  Ep, etc.
%
%  [Ep,Z] = Zpredict(Z,B,fe) interpolates Z using 
%  Zinterp(fe,Z,fftfreqp(size(B,1)) prior to computing Ep.
%
%  [Ep,Z] = Zpredict(Z,B,fe,opts)
%
%  See also ZPREDICT_DEMO.

verbose = 0;

N  = size(B,1);
[~,fg] = fftfreq(N); % Will include f=0.5 if N is even.

if nargin > 2

    if verbose
      fprintf('Zpredict.m: Interpolating Z onto same frequency grid as input.\n');
    end
    
    if nargin < 4
        Z = Zinterp(fe,Z,fg);
    else
        Z = Zinterp(fe,Z,fg,opts);
    end

    if mod(N,2) == 0
      % Last element of Z from Zinterp corresponds to fg = 0.5 if N is even.
      Z = [ Z(1:end,:) ; flipud(conj(Z(2:end-1,:))) ];
    else
      Z = [ Z(1:end,:) ; flipud(conj(Z(2:end,:))) ];
    end

end

assert(size(B,1) == size(Z,1),'Number of rows in B must equal number of rows in Z');

Nout = size(Z,2)/size(B,2); % Number of output columns
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

% TODO: Verify that complex part is small. (May be non-zero due to
% roundoff?)

Ep = real(Ep);
