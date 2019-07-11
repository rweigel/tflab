function opts = transferfnFD_options(rn)

opts = struct();

opts.filestr = sprintf('transferfnFD_options-%d',rn);

opts.info = struct();
    opts.info.instr = 'B'; % Or {'$B_x$',...} (1 cell element per column in In)
    opts.info.inunit= 'nT';
    opts.info.outstr = 'E'; % Or {'$E_x$',...} (1 cell element per column in Out)
    opts.info.outunit= 'mV/m';
    opts.info.timeunit = 's';
    opts.info.timestr = 'Time since 2019-01-01';

opts.transferfnFD = struct();
    opts.transferfnFD.log = 0;
    opts.transferfnFD.plot = struct();
        opts.transferfnFD.plot.timeseries = [0,0,0];
        opts.transferfnFD.plot.spectrum = [0,0,0];
        opts.transferfnFD.plot.Z = [0,0,0];
        opts.transferfnFD.plot.H = [0,0,0];

% # of points at start and end to trim before computing
% pe/cc/mse/sn/coherence
opts.td.Ntrim = NaN;

opts.td.detrend.function = @removemean;
opts.td.detrend.functionstr = '';
opts.td.detrend.functionargs = {};

opts.td.dt    = 1;  % Dimensionless time step
opts.td.start = 1;  % Dimensionless start time
% opts.td.start can also be a time string of the form
% 'yyyy-mm-ddTHH:MM:SS.FFF'. In this case, opts.td.dt units are
% milliseconds. Example: 
% opts.td.start = '2001-01-01T00:00:00.000';
% opts.td.dt    = 1000; % 1-second cadence.

opts.td.window = struct();
    opts.td.window.function = ''; 
    opts.td.window.functionstr = 'none';
    opts.td.window.functionargs = {};
    opts.td.window.width = NaN;
    opts.td.window.shift = NaN;
    opts.td.window.plot = [0,0,0];
    opts.td.window.log = 0;

if 0
    % Note: Same window applied to input and output
    % See 'help rectwin' for list of available functions.
    opts.td.window.function = @tdwindow; 
    opts.td.window.functionstr = 'Rectangular';
    % functionargs are arguments after first argument to function, which is the
    % length of the segment.
    opts.td.window.functionargs = {@rectwin};
end

opts.td.prewhiten = struct();
    opts.td.prewhiten.function = '';
    opts.td.prewhiten.plot = [0,0,0];
    opts.td.prewhiten.log = 0;

if 0
    % Note: Same prewhitening filter applied to input and output
    opts.td.prewhiten.function = @tdprewhiten;
    opts.td.prewhiten.functionstr = 'diff';
    opts.td.prewhiten.functionargs = {};
end

opts.fd.evalfreq = struct();
    opts.fd.evalfreq.plot = [0,0,0];
    opts.fd.evalfreq.log = 0;

    opts.fd.evalfreq.function = @evalfreq;   
    % The following args are the arguments to evalfreq after the first
    % argument which is the number of time points (and is computed
    % internally).
    opts.fd.evalfreq.functionargs = {7, 'logarithmic'};
    opts.fd.evalfreq.functionstr  = ...
        sprintf('%d frequencies per decade',...
                opts.fd.evalfreq.functionargs{1});

opts.fd.window = struct();
    opts.fd.window.function = @rectwin; 
    opts.fd.window.functionstr = 'rectangular';
    opts.fd.window.log = 0; 
    opts.fd.window.plot = 0; 

opts.fd.stack = struct();
    opts.fd.stack.average = struct();
        opts.fd.stack.average.function = @transferfnAverage;
        opts.fd.stack.average.functionstr = '';
        opts.fd.stack.average.functionargs = {};

opts.fd.interpolation = struct();
    opts.fd.interpolation.function = @interp1;
    opts.fd.interpolation.functionstr = 'interp1()';
    opts.fd.interpolation.functionargs = {'linear',0};
    %opts.fd.interpolation.functionargs = {'linear','extrap'};
    
opts.fd.regression = struct();
    %opts.fd.regression.function = @ols_analytic;
    %opts.fd.regression.functionstr = 'OLS using analytic formula';
    %opts.fd.regression.functionargs = {};

    opts.fd.regression.function = @ols_regress;
    opts.fd.regression.functionstr = 'OLS using regress() function';
    opts.fd.regression.functionargs = {};
    opts.fd.regression.plot = 0;
    opts.fd.regression.log = 0;

    %opts.fd.regression.function = @robust_robustfit;
    %opts.fd.regression.functionstr = 'Robust regression using robustfit() function';
    %opts.fd.regression.functionargs = {[],[],'off'};

    %opts.fd.regression.function = @robust_v1;
    %opts.fd.regression.functionstr = 'Robust regression using robust_v1() function';
    %ropts = struct();
    %ropts.weightfn = 'huber';
    %ropts.stepmax = 50;
    %ropts.zeps = sqrt(eps);
    %ropts.hardcut = Inf; % 2.8
    %ropts.snstop = 1000;
    %ropts.verbose = 0;    
    %opts.fd.regression.functionargs = {ropts};

if rn == 0
    % When no noise, should get exact TF used to generate
    % the output data (within limits of numerical precision).
    opts.description = 'Default except dN = 0, w = 0.';
    opts.fd.evalfreq.functionargs = {[1,0], 'linear'};
    opts.fd.evalfreq.functionstr  = 'dN = 0, w = 0';
elseif rn == 1
    opts.description = 'Default';
elseif rn == 2
    opts.description = 'yulewalker(10) prewhiten';
    opts.td.prewhiten.method = 'yulewalker';
    opts.td.prewhiten.methodstr = 'yulewalker10';
    opts.td.prewhiten.options = 10;
elseif rn == 3
    opts.description = 'Parzen window in FD';
    opts.fd.window.function = @parzenwin; 
    opts.fd.window.functionstr = 'parzen';
elseif rn == 4
    opts.description = 'Robust regression';
    opts.fd.regression.method = 'robust_robustfit'; 
    opts.fd.regression.methodstr = 'Robust regression using robustfit()';
elseif rn == 5
    opts.description = 'PCA rotation';
    opts.td.transform = 'pca';
    %opts.td.transform.methodstr = 'pca';
elseif rn == 6
    opts.description = '1.5 day window';
    opts.td.window.width = 3600*48;
    opts.td.window.shift = 3600*24;
else
    error('Invalid option set number');
end

