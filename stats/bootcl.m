function lims = bootcl(z,fn,N,f)
%BOOTCL Bootstrap confidence limits
%
%   lims = BOOTCL(z) - for each row in z, compute N=1000 bootstrap means
%   by resampling 63% of the columns with replacement. Each row of lims
%   contains the value nearest the lowest/highest 2.5%/97.5% of the
%   N boostrap means.
%
%   lims = BOOTCL(z,fn)
%   lims = BOOTCL(z,fn,N)
%   lims = BOOTCL(z,fn,N,f)
%

% TODO: Verify that 63% is used.

z = squeeze(z);

if nargin < 2
    fn = @mean;
end
if nargin < 3
    N = 1000;
end
if nargin < 4
    f = 0.025;
end

n = round(N*f);

z = transpose(z);

% Each row in V contains N bootstrap samples from that row in z
V = bootstrp(N,fn,z); 
V = sort(V,1); % Sort each column
l = V(n,:);     % Select the nth lowest value
u = V(N-n+1,:); % Select the N-n th highest value
lims = [l',u'];

end
