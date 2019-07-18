function [h,t,zg,fg] = Z2H(Z,fe,N,verbose)
% Z2H Convert frequency domain transfer function to impulse response
%
%   [h,t] = Z2H(Z) returns h = ifft(Z) and t = N*fftfreq(N) where N is
%   number of elements in first non-singleton dimension.
% 
%   Z2H(Z,f) uses f to compute Z for negative frequencies
% 
%   Z2H(Z,f,N) interpolates Z onto N-point DFT grid using Zi =
%   ZINTERP(Z,f,fftfreqp(N)) prior to computing h. f must contain only
%   postive DFT frequencies and may contain f = 0.5.
%
%   Z2H(Z,f,fg) interpolates Z onto frequency grid given by array fg using
%   Zi = ZINTERP(f,Z,fg) prior to computing h. f and fg must contain only
%   postive DFT frequencies and may contain f=0.5.
%
%   Z2H(f,Z,N,verbose)
%
%   See also H2Z, ZINTERP.


if nargin < 4
    verbose = 0;
else
    addpath([fileparts(mfilename('fullpath')),'/logging']);
end

verbose = 1;

assert(nargin ~= 2,'Number of inputs must be 1, 3, or 4.');
assert(ndims(Z) <= 2,'Required: ndims(Z) <= 2');

if nargin == 1
    h = ifft(Z);
    if numel(Z) == length(Z) && size(Z,2) > 1
        N = size(Z,2);
        t = N*fftfreq(N);
    else
        N = size(Z,1);
        t = round(N*fftfreq(N))';
    end
    return;
end

flip = 0;
if numel(Z) == length(Z) && size(Z,2) > 1
    assert(length(Z) == length(fe),...
        'length(Z) == length(fe) is required when ndims(Z) == 1.');
    Z = transpose(Z);
    flip = 1;
end
assert(size(Z,1) == length(fe),'size(Z,1) == length(fe) is required');

if nargin > 2
    if verbose
        logmsg(dbstack,...
            'Interpolating Z onto frequency grid with spacing of 1/%d\n',...
            N);
    end
    
    % fg = Unique fft frequencies (with -0.5 mapped to +0.5)
    [~,fg] = fftfreq(N);
    Zi = Zinterp(fe,Z,fg');

    if verbose == 2
        for i = 1:size(Z,2)
            figure()
            loglog(fe,real(Z(:,1)));
            hold on;grid on;box on;
            loglog(fe,imag(Z(:,1)));
        end
    end
end

if mod(N,2) == 0
    % Last fg value is f=0.5. This moves Z row into proper place.
    % No conj. needed because it is real.
    Zfull = [Zi(1:end,:) ; flipud(conj(Zi(2:end-1,:)))];
else
    Zfull = [Zi(1:end,:) ; flipud(conj(Zi(2:end,:)))];
end

zg = Zfull;
fg = fftfreq(size(Zfull,1))';

% Compute impulse response
h = ifft(Zfull);
t = N*fftfreq(N)';

assert(isreal(h),...
            ['Computed impulse response has imaginary component.'...
            'Check calculation of Z']);

if flip
    h = h';
    t = t';
    zg = transpose(zg);
    fg = fg';
end