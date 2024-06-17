diary bootcl_demo.log
n = 100;
x = randn(1,n);
lims = bootcl(x);

fprintf('Bootstrap: [%.4f, %.4f]\n', lims);
fprintf('Exact:     [%.4f, %.4f]\n',-1.96/sqrt(n),1.96/sqrt(n));

% Test with 10 rows. Bootstrap performed on each row.
x = randn(10,n);
lims = bootcl(x);
fprintf('\nBootstrap: [%.4f, %.4f]', lims');
fprintf('\nAverages:  [%.4f, %.4f]\n', mean(lims));
diary off