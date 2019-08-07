function [Zi,fi] = zinterp(f,Z,fi,varargin)
% ZINTERP - Interpolate transfer function onto frequency grid
%
%  [Zi,fi] = zinterp(f,Z,N) returns Zi on N-point DFT grid fi = fftfreq(N)
%  using INTERP1. f >= 0 is required and f(end) = +0.5 is allowed.
%
%  If f(1) = 0, it is not used for interpolation and Zi(1) is set to Z(1)
%  (or Zi(1,:) = Z(1,:) if Z is a matrix).
%
%  Zi = zinterp(f,Z,fi) returns interpolated values of Zi at frequencies fi
%  using the given Z at frequencies f using INTERP1. fi >= 0 and f >= 0 are
%  required.
%
%  If f(1) = 0 and fi(1) = 0, Zi(1) is set to Z(1) (or Zi(1,:) = Z(1,:) if
%  Z is a matrix). the first row of Zg is zeros.
%
%  If f(end) = +0.5, the interpolation for the imaginary part of Z is
%  performed using only Z(1:end-1,:) and imag(Zg(end,:)) is set to zeros.
%
%  See also ZINTERP_DEMO, ZINTERP_TEST.

% TODO: Allow interpolation in log space.

addpath([fileparts(mfilename('fullpath')),'/../misc/']);
verbose = 1;

if nargin < 4
    varargin = {'linear',0};
else
    if isstruct(varargin{1})
        varargin = varargin{1}.fd.interpolation.functionargs;
    end
end

assert(ndims(Z) == 2,'Z can have at most two dimensions.');

transposeZ = 0;
if size(Z,1) == 1
    assert(size(Z,2) == length(f),'If Z is a vector, it must have same length as f.');
    transposeZ = 1;
    Z = transpose(Z);
else
    assert(size(Z,1) == length(f),'If Z is a matrix, number of rows must equal length(f).');
end

N = NaN;
if isscalar(fi)
    assert(fi > 1,'N > 1 is required when using zinterp(f,Z,N)');
    % Zinterp(f,Z,N) usage
    N = fi;
    [fa,fi] = fftfreq(N);
    % fi contains unique DFT frequencies with f = -0.5 mapped to f = +0.5
    % when N is even so that fg(end) = +0.5.
end

if any(f < 0) || any(fi < 0)
    error('Elements of f and fi must be greater than or equal to zero.');
end

if nargin == 3 && (length(f) == length(fi))
    if all(f(:) == fi(:))
        if verbose
            logmsg(dbstack,...
                ['all(f == fi) returned true. '...
                 'No interpolation will be performed.\n']);
        end
        if ~isnan(N)
            % Zinterp(f,Z,N) usage.
            [Zi,fi] = Zfull(Z,N);
        else
            Zi = Z;
        end
        return;
    end
end

if verbose
    logmsg(dbstack, 'First interp frequency: %.4f\n',fi(1));
    logmsg(dbstack, 'First given frequency : %.4f\n',f(1));
    logmsg(dbstack, 'Last interp frequency : %.4f\n',fi(end));
    logmsg(dbstack, 'Last given frequency  : %.4f\n',f(end));
end

% Remove Z value for fe = 0 if found.
fe0 = 0;
if f(1) == 0
    fe0 = 1;
    f = f(2:end);
    Zo = Z(1,:);
    Z = Z(2:end,:);
    assert(all(imag(Zo)) == 0,'If f(1) == 0 expect all(imag(Z(1,:))==0)');
end

% If fg(1) == 0, remove it before interpolation.
fg0 = 0; 
if fi(1) == 0
    fg0 = 1;
    fi = fi(2:end);
end

if f(end) == 0.5
    assert(all(imag(Z(end,:)) == 0),'If f(end) == 0.5 expect all(imag(Z(end,:))==0)');
end

Zi = zeros(length(fi),size(Z,2));

for k = 1:size(Z,2)
    
    Ig = ~isnan(Z(:,k)); % "good" values
    if verbose && (length(Ig) ~= size(Z,1))
        logmsg(dbstack, ...
            'Dropping %d NaN values in Z(:,%d)\n',...
            size(Z,1)-length(Ig),k);
    end
    
    fek = f(Ig);
    Zk = Z(Ig,k);

    if f(end) ~= 0.5
        Zi(:,k) = interp1(fek,Zk,fi,varargin{:});
    else
        % Z at f = 0.5 should have no imaginary component.
        % This performs imaginary interpolation only for 0 < f < 0.5.
        Zir = interp1(fek,real(Zk),fi,varargin{:});
        Zii = interp1(fek(1:end-1),imag(Zk(1:end-1)),fi,varargin{:});
        Zi(:,k) = Zir + sqrt(-1)*Zii;
    end

end

if fi(end) == 0.5
    Zi(end,:) = real(Zi(end,:));
end

if fg0 % If lowest grid frequency is zero (and so was removed)
    if fe0
        % If Z was given at fe = 0, use it as the "interpolated" value.
        Zi = [Zo;Zi];
    else
        % Otherwise, set it to zero
        Zi = [zeros(1,size(Z,2));Zi];
    end
end

if ~isnan(N)
    % Zinterp(f,Z,N) usage.
    [Zi,fi] = Zfull(Zi,N);
end

if transposeZ
    Zi = transpose(Zi);
end

function [Zi,fi] = Zfull(Zi,N)
    if mod(N,2) == 0
        % Last Zi row corresponds to fg = +0.5. This put Z in conventional
        % frequency order for even N where
        % f = [0, 1/N, ..., N/2-1 , -0.5, ..., -1/N]
        Zi = [Zi(1:end,:) ; flipud(conj(Zi(2:end-1,:)))];
    else
        Zi = [Zi(1:end,:) ; flipud(conj(Zi(2:end,:)))];
    end
    fi = fftfreq(N)';
end

end
