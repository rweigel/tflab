function [bands, f, fe] = dftbands(x, Ic, Ne)
%DFTBANDS Compute DFT and split into bands
%
%   [bands, f, fe] = DFTBANDS(x, Ic, Ne) computes xdft = fftu(x)
%   and splits the result into bands determined by
%   bands{j,1} = xdft(Ic(j)-Ne(j):Ic(j)+Ne(j),:)
%
%   Example:
%
%     x = randn(100,2);
%     [fc, Ic, Ne] = evalfreq(size(x,1));
%     [bands, f] = dftbands(x, Ic, Ne)
%     whos bands f
%
%   [bands, f, fc] = DFTBANDS(x, opts) uses opts.fd.evalfreq.function with
%   arguments opts.fd.evalfreq.functionargs to compute Ic and Ne.
%    
%   See also evalfreq.

if nargin == 2
    opts = Ic;
    [fe,Ic,Ne] = opts.fd.evalfreq.function(...
                        size(x,1),opts.fd.evalfreq.functionargs{:});
end

[dftu, fu] = fftu(x);

for j = 1:length(Ic)
    r = Ic(j)-Ne(j):Ic(j)+Ne(j); % Index range
    bands{j,1} = dftu(r,:);
    f{j,1} = fu(r,:);
end
