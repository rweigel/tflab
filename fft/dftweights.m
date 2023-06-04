function W = dftweights(f, dftE, dftB, opts)

winfn = opts.fd.window.function;
if opts.fd.window.loglevel
    logmsg('Using FD window function: %s\n',func2str(winfn));
end

for j = 1:length(f)
    w = winfn(length(f{j}));
    W{j,1} = w/sum(w);
end


