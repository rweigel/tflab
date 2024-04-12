function [wE,wB] = dftweights(f, dftIn, dftOut, opts)
%DFTWEIGHTS  Compute weights for DFT coefficients
%
%  [wIn,wOut] = DFTWEIGHTS(f, dftIn, dftOut, opts) computes weights using
%  opts.fd.window.function. f, dftIn, and dftOut are cell arrays.
%
%  See also TDWINDOW.

winfn = opts.fd.window.function;

for s = 1:length(f)
    w = winfn(length(f{s}));
    for j = 1:size(dftIn{s},2)
        wE{s,1}(:,j) = w/sum(w);
    end
    for j = 1:size(dftOut{s},2)
        wB{s,1}(:,j) = w/sum(w);
    end
end