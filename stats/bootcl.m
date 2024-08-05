function lims = bootcl(x,fn,N,p,use_bootstrp)
%BOOTCL Bootstrap confidence limits
%
%   LIMS = BOOTCL(X) compute 1000 bootstrap means for each column of X by
%   resampling it with replacement.
%
%   LIMS(1,:) contains the 2.3% bootstrap percentile value (lower limit)
%   LIMS(2,:) contains the 97.7% bootstrap percentile value (upper limit)
%
%   PCTILE(...,'Method','exact') is used to compute the percentiles
%
%   lims = BOOTCL(X,FN)
%   lims = BOOTCL(X,FN,N)
%   lims = BOOTCL(X,FN,N,F)
%

x = squeeze(x);

s = size(x);
assert(length(s) <= 2,'length(size(squeeze(X)) must be <= 2')

flip = 0;
if s(1) == 1
    flip = 1;
    x = transpose(x);
end

if nargin < 2
    fn = @mean;
end
if nargin < 3
    N = 1000;
end
if nargin < 4
    p = normcdf(-2);
    p = 100*[p, 1-p];
end
if nargin < 5
    use_bootstrp = 1;
end

% Each row in xb contains N bootstrap samples from that row in x
if ~isempty(ver('stats')) && use_bootstrp == 1
    xb = bootstrp(N,fn,x);
else
    xb = zeros(N, size(x,2));
    nr = size(x,1);
    for i = 1:size(x,2)
        for b = 1:N
            I = randsample(nr,nr,1);
            xb(b,i) = fn(x(I,i));
        end
    end
end

l = prctile(xb,p(1));
u = prctile(xb,p(2));

% TODO: Use pivot
% https://www2.stat.duke.edu/~banks/111-lectures.dir/lect13.pdf

if flip == 0
    lims = [l; u];
else
    lims = [l, u];
end