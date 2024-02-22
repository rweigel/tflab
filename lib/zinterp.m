function [Zi,fi,Zip,fip] = zinterp(f,Z,fi,interp1args)
% ZINTERP - Interpolate transfer function on to frequency grid
%
%  [Zi,fi] = ZINTERP(f,Z,N) returns Zi on the N-point DFT grid fi given by
%  fi = fftfreq(N) using INTERP1 with method='linear' and extrapval=0. f >=
%  0 is required and f(end) = +0.5 is allowed. The ordering of values in fi
%  follows the convention of FFTFREQ (and FFT) and a given Z value at f =
%  +0.5 will appear at fi = -0.5 if N is even.
%
%  [Zi,fi,Zip,fip] = ZINTERP(f,Z,N) returns Zip for positive frequencies
%  fip.
%
%  If f(1) = 0, it is not used for interpolation and if fi(1) = 0, Zi(1,:)
%  is set to Z(1,:).
%
%  Zi = ZINTERP(f,Z,fi) returns interpolated values of Zi at frequencies fi
%  using the given Z at frequencies f using INTERP1. fi >= 0 and f >= 0 are
%  required and the f(1) = 0 case is handled as described above.
%
%  See also ZINTERP_DEMO, ZINTERP_TEST.

% TODO: Allow interpolation in log space.

if nargin < 4
    interp1args = {'linear',0};
end

assert(ismatrix(Z),'Z can have at most two dimensions.');
assert(iscolumn(f),'f must be a column vector (nx1)');
assert(iscolumn(fi),'fi must be a column vector (nx1)');
assert(size(Z,1) == size(f,1),'Z and f must have same number of rows');

% Keep any row of Z without a NaN
% TODO: Handle cases where too few values left for interpolation.
Ik = any(~isnan(Z),2);
f = f(Ik);
Z = Z(Ik,:);

N = NaN;
if isscalar(fi) % Zinterp(f,Z,N) usage
    assert(fi > 1,'N > 1 is required when using zinterp(f,Z,N)');
    N = fi;
    [fa,fi] = fftfreq(N);
    % fi contains unique DFT frequencies with f = -0.5 mapped to f = +0.5
    % when N is even so that fi(end) = +0.5.
end

if length(f) == 1
    I = find(f == fi,1);
    if ~isempty(I)
        logmsg('Only one frequency that matches an interpolation frequency.\n');
        logmsg('Setting Z = 0 at all interpolation frequencies execept matching frequency.\n');
        Zo = Z;
        Z = zeros(length(fi),size(Zo,2));
        Z(I,:) = Zo;
        f = fi;
    else
        error('Only one frequency and it does not equal any interpolation frequency. Cannot interpolate.');
    end
end

if any(f < 0) || any(fi < 0)
    error('Elements of f and fi must be greater than or equal to zero.');
end

if length(f) == length(fi) && all(f(:) == fi(:))
    logmsg('all(f == fi) returned true. No interpolation will be performed.\n');
    if ~isnan(N)
        % Zinterp(f,Z,N) usage.
        fir = fi;
        Zir = Z;
        % Create Z having negative frequency elements.
        [Zi,fi] = zfull(Z,N);
    else
        Zi = Z;
    end
    return;
end

if 0
    logmsg('First interp frequency: %.4f\n',fi(1));
    logmsg('First given frequency : %.4f\n',f(1));
    logmsg('Last interp frequency : %.4f\n',fi(end));
    logmsg('Last given frequency  : %.4f\n',f(end));
end

% Remove Z value for f = 0 if found.
f0 = 0;
if f(1) == 0
    f0 = 1;
    f = f(2:end);
    Z0 = Z(1,:);
    Z = Z(2:end,:);
    assert(all(imag(Z0)) == 0,'If f(1) == 0 expect all(imag(Z(1,:))==0)');
end

if f(end) == 0.5
    % TODO: This assumes f is normalized
    assert(all(imag(Z(end,:)) == 0),'If f(end) == 0.5 expect all(imag(Z(end,:))==0)');
    Zi_re = interp1(f,real(Z),fi,interp1args{:});
    Zi_im = interp1(f(1:end-1),imag(Z(1:end-1,:)),fi,interp1args{:});
    Zi = Zi_re + 1j*Zi_im;
else
    Zi = interp1(f,Z,fi,interp1args{:});
end

if fi(1) == 0 % If lowest interp. frequency is zero (and so was removed)
    if f0 
        % If Z was given at f = 0, use it as the "interpolated" value.
        Zi(1,:) = Z0;
    else
        % Otherwise, set it to extrapval used in interp1 call.
        if isscalar(interp1args{2})
            % If 'extrap' is given instead of an extrap value, e.g.,
            % {'linear','extrap'}.
            Zi(1,:) = interp1args{2}*ones(1,size(Z,2));
        end
    end
end

if ~isnan(N)
    Zip = Zi; % Zi for positive frequencies only
    fip = fi;
    % Zinterp(f,Z,N) usage. Create z with negative frequencies (for use
    % with ifft, for example).
    [Zi,fi] = zfull(Zi,N);
end

function [Zi,fi] = zfull(Zi,N)
    if mod(N,2) == 0
        % Last Zi row corresponds to fi = +0.5 (or equivalent un-normalized). 
        % This puts Z in conventional frequency order for even N where
        % f = [0, 1/N, ..., N/2-1 , -0.5, ..., -1/N]
        Zi = [Zi(1:end,:) ; flipud(conj(Zi(2:end-1,:)))];
        % Force f = +0.5 (or equivalent un-normalized) element to be real  
        Zi(end/2+1,:) = real(Zi(end/2+1,:));
    else
        Zi = [Zi(1:end,:) ; flipud(conj(Zi(2:end,:)))];
    end
    fi = fftfreq(N)';
end

end
