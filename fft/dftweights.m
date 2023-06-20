function [wE,wB] = dftweights(f, dftE, dftB, opts)

winfn = opts.fd.window.function;

for s = 1:length(f)
    w = winfn(length(f{s}));
    for j = 1:size(dftE{s},2)
        wE{s,1}(:,j) = w/sum(w);
    end
    for j = 1:size(dftB{s},2)
        wB{s,1}(:,j) = w/sum(w);
    end
end


