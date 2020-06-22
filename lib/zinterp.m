function [Zi,fi] = zinterp(f,Z,fi,opts)
% ZINTERP - Interpolate transfer function onto frequency grid
%
%  [Zi,fi] = ZINTERP(f,Z,N) returns Zi on N-point DFT grid fi = fftfreq(N)
%  using INTERP1 with method='linear' and extrapval=0. f >= 0 is required
%  and f(end) = +0.5 is allowed. The ordering of values in fi follows the
%  convention of FFTFREQ (and FFT) and a given Z value at f = +0.5 will
%  appear at fi = -0.5 if N is even.
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

addpath([fileparts(mfilename('fullpath')),'/../misc/']);

if nargin < 4
    opts = struct('loglevel',0,'interp1args',{{'linear',0}});
else
    assert(isstruct(opts),'opts must be a structure');
    if ~isfield(opts,'loglevel')
        opts.loglevel = 0;
    end
    if ~isfield(opts,'interp1args')
        opts.interp1args = {'linear',0};
    end
end

assert(ndims(Z) == 2,'Z can have at most two dimensions.');

assert(iscolumn(f),'f must be a column vector (nx1)');
assert(iscolumn(fi),'fi must be a column vector (nx1)');
assert(size(Z,1) == size(f,1),'Z and f must have same number of rows');

N = NaN;
if isscalar(fi) % Zinterp(f,Z,N) usage
    assert(fi > 1,'N > 1 is required when using zinterp(f,Z,N)');
    N = fi;
    [fa,fi] = fftfreq(N);
    fa = fa';
    fi = fi';
    % fi contains unique DFT frequencies with f = -0.5 mapped to f = +0.5
    % when N is even so that fi(end) = +0.5.
end

if any(f < 0) || any(fi < 0)
    error('Elements of f and fi must be greater than or equal to zero.');
end

if length(f) == length(fi) && all(f(:) == fi(:))
    if opts.loglevel > 0
        logmsg(...
            ['all(f == fi) returned true. '...
             'No interpolation will be performed.\n']);
    end
    if ~isnan(N)
        % Zinterp(f,Z,N) usage.
        % Create Z having negative frequency elements.
        [Zi,fi] = zfull(Z,N);
    else
        Zi = Z;
    end
    return;
end

if opts.loglevel > 0
    logmsg( 'First interp frequency: %.4f\n',fi(1));
    logmsg( 'First given frequency : %.4f\n',f(1));
    logmsg( 'Last interp frequency : %.4f\n',fi(end));
    logmsg( 'Last given frequency  : %.4f\n',f(end));
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
    % Can't do this because we don't know if f is normalized and generally
    % don't know normalization.
    %assert(all(imag(Z(end,:)) == 0),'If f(end) == 0.5 expect all(imag(Z(end,:))==0)');
end

Zi = interp1(f,Z,fi,opts.interp1args{:});

if fi(1) == 0 % If lowest interp. frequency is zero (and so was removed)
    if f0 
        % If Z was given at f = 0, use it as the "interpolated" value.
        Zi(1,:) = Z0;
    else
        % Otherwise, set it to extrapval used in interp1 call.
        if isscalar(opts.interp1args{2})
            % If 'extrap' is given instead of an extrap value, e.g.,
            % {'linear','extrap'}.
            Zi(1,:) = opts.interp1args{2}*ones(1,size(Z,2));
        end
    end
end

if ~isnan(N)
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
