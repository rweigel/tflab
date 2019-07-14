function Zi = Zinterp(fe,Z,fg,varargin)
% ZINTERP - Interpolate transfer function onto frequency grid
%
%  Zg = Zinterp(fe,Z,fg) returns Zg containing interpolated values of Z at
%  frequencies fg, where fe >= 0 and fg >= 0 are evaluation and grid
%  frequencies, respectively. Rows of Z that have NaNs are ignored. fg
%  values that are outside of range of fe values have an interpolated Zg
%  value of zero (instead of NaN as would be the case for interp1).
%
%  Elements of fe are frequencies of corresponding elements of Z if Z is a
%  vector. If Z is a matrix, each column of Z is interpolated.
%
%  If fe(1) = 0, it is not used for interpolation. If fg(1) = 0 and fe(1) =
%  0, the first row of Zg is set equal to first row of Z; otherwise, the
%  first row of Zg is zeros.
%
%  If fe(end) = 0.5, the interpolation for the imaginary part of Z is
%  performed using only Z(1:end-1,:) and imag(Zg(end,:)) is set to zeros.
%
%  Zg = Zinterp(fe,Z,fg,opts)
%
%  Zg = Zinterp(fe,Z,fg,)
%
%  See also ZINTERP_DEMO.

% TODO: Allow interpolation in log space.

verbose = 0;

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
    assert(size(Z,2) == length(fe),'If Z is a vector, it must have same length as fe.');
    transposeZ = 1;
    Z = transpose(Z);
else
    assert(size(Z,1) == length(fe),'If Z is a matrix, number of rows must equal length(fe).');
end

if any(fe < 0) || any(fg < 0)
    error('Elements of fe and fg must be greater than or equal to zero.');
end

if nargin == 3 && (length(fe) == length(fg))
    if all(fe(:) == fg(:))
        if verbose
            warning('all(fe == fg) returned true.  No interpolation will be performed.');
        end
        Zi = Z;
        return;
    end
end

if verbose
    fprintf('Zinterp.m: First grid frequency       : %.4f\n',fg(1));
    fprintf('Zinterp.m: First evaluation frequency : %.4f\n',fe(1));
    fprintf('Zinterp.m: Last grid frequency        : %.4f\n',fg(end));
    fprintf('Zinterp.m: Last evaluation frequency  : %.4f\n',fe(end));
end

% Remove fe = 0
fe0 = 0; % Will be 1 if f=0 fe value is zero.
if fe(1) == 0
    fe0 = 1;
    fe = fe(2:end);
    Zo = Z(1,:);
    Z = Z(2:end,:);
    assert(all(imag(Zo)) == 0,'If fe(1) == 0 expect all(imag(Z(1,:))==0)');
end

% Remove fg = 0
fg0 = 0; % Will be 1 if f = 0 fg value is zero.
if fg(1) == 0
    fg0 = 1;
    fg = fg(2:end);
end

if fe(end) == 0.5
    assert(all(imag(Z(end,:)) == 0),'If fe(end) == 0.5 expect all(imag(Z(end,:))==0)');
end

Zi = zeros(length(fg),size(Z,2));

for k = 1:size(Z,2)
    
    Ig = ~isnan(Z(:,k)); % "good" values
    if verbose && (length(Ig) ~= size(Z,1))
        fprintf('Zinterp.m: Dropping %d NaN values in Z(:,%d)\n',size(Z,1)-length(Ig),k);
    end
    
    fek = fe(Ig);
    Zk = Z(Ig,k);

    if fe(end) ~= 0.5
        Zi(:,k) = interp1(fek,Zk,fg,varargin{:});
    else
        % Z at f = 0.5 should have no imaginary component.
        % This performs imaginary interpolation only for 0 < f < 0.5.
        Zir = interp1(fek,real(Zk),fg,varargin{:});
        Zii = interp1(fek(1:end-1),imag(Zk(1:end-1)),fg,varargin{:});
        Zi(:,k) = Zir + sqrt(-1)*Zii;
    end

end

if fg(end) == 0.5
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

if transposeZ
    Zi = transpose(Zi);
end
