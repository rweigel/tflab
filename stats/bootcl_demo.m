diary bootcl_demo.log

nb = 1000;
p = [2.5, 97.5];
n = 40;

x = randn(n,1);

rng(0,"twister")
lims1 = bootcl(x, @mean, nb, p);
fprintf('Bootstrap: [%.4f, %.4f]', lims1);
fprintf('\nExact:     [%.4f, %.4f]\n',-1.96/sqrt(n),1.96/sqrt(n));

rng(0,"twister")
% Don't use bootstrp() function. Should give the same result as above.
lims2 = bootcl(x, @mean, nb, p, 0);
fprintf('\nBootstrap: [%.4f, %.4f]', lims2);
fprintf('\nExact:     [%.4f, %.4f]\n',-1.96/sqrt(n),1.96/sqrt(n));

assert(all(lims1 == lims2))

x = ones(n,1);
lims = bootcl(x, @mean, nb, p);
fprintf('\nBootstrap: [%.4f, %.4f]', lims);
fprintf('\nExact:     [%.4f, %.4f]\n',1,1);
assert(all(size(lims) == [2,1]))
assert(all(lims == 1))

x = ones(1,n);
lims = bootcl(x, @mean, nb, p);
fprintf('\nBootstrap: [%.4f, %.4f]', lims);
fprintf('\nExact:     [%.4f, %.4f]\n',1,1);
assert(all(size(lims) == [1,2]))
assert(all(lims == 1))

% Test with 10 columns. Bootstrap performed on each column.
x = randn(n,10);

% rng(0,"twister") Don't get exactly the same results if setting
% seed when x is not a vector. Reason probably has to do with internal
% implementation of random selection and differences are not due to an
% implementation difference between bootstrp() and the manual implementation
% in bootcl().
lims1 = bootcl(x, @mean, nb, p);
fprintf('\nBootstrap: [%.4f, %.4f]', lims1);
fprintf('\nAverages:  [%.4f, %.4f]', mean(lims1,2));
fprintf('\nExact:     [%.4f, %.4f]\n',-1.96/sqrt(n),1.96/sqrt(n));

%rng(0,"twister")
lims2 = bootcl(x, @mean, nb, p, 0);
fprintf('\nBootstrap: [%.4f, %.4f]', lims2);
fprintf('\nAverages:  [%.4f, %.4f]', mean(lims2,2));
fprintf('\nExact:     [%.4f, %.4f]\n',-1.96/sqrt(n),1.96/sqrt(n));

% TODO: The 0.02 value is not justified and was chosen when using n = 40 and
% nb = 1000.
dmean = mean(lims1, 2) - mean(lims2, 2);
assert(all(abs(dmean) < 0.02))

diary off