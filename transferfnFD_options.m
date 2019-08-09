function opts = transferfnFD_options(os)
%TRANSFERFNFD_OPTIONS - Return options for transferfnFD().
%
%  opts = TRANSFERFNFD_OPTIONS() returns default options.
%
%  opts = TRANSFERFNFD_OPTIONS(os), where os is an option set.
%
%  See also TRANSFERFNFD, TRANSFERFNFD_DEMO.

opts = struct();

opts.filestr = sprintf('transferfnFD_options-%d',os);

opts.info = struct();
    opts.info.instr = 'B'; % Cell array or string. 
    % Or {'$B_x$', ...} (1 cell element per column in In)
    
    opts.info.outstr = 'E'; % Cell array or string. 
    % Or {'$E_x$', ...} (1 cell element per column in Out)
    
    % Used if time array not passed to transferfnFD.
    opts.info.timestr = 'Time since 2000-01-01'; 

    opts.info.inunit= 'nT';
    opts.info.outunit= 'mV/m';
    opts.info.timeunit = 's';
    
opts.transferfnFD = struct();
    opts.transferfnFD.loglevel = 1;
    % Elements of arrays for plot.* are 0 or 1 and
    % [showplot, savepng, savepdf] 
    opts.transferfnFD.plot = struct();
        opts.transferfnFD.plot.timeseries = [0,0,0];
        opts.transferfnFD.plot.spectrum = [0,0,0];
        opts.transferfnFD.plot.Z = [0,0,0];
        opts.transferfnFD.plot.H = [0,0,0];

% # of points at start and end to trim before computing metrics
% (pe/cc/mse/sn/coherence)
opts.td.Ntrim = NaN;

opts.td.detrend.function = @removemean;
opts.td.detrend.functionstr = '';  % Optional descriptive name
opts.td.detrend.functionargs = {}; % Arguments after first argument to fn.

% Dimensionless time; ignored if time array passed to transferfnFD.
opts.td.dt    = 1;  
% Dimensionless start; ignored if time array passed to transferfnFD.
opts.td.start = 1; 

% opts.td.start can also be a time string of the form
% 'yyyy-mm-ddTHH:MM:SS.FFF'. In this case, opts.td.dt units are milliseconds.
% Example: 
%   opts.td.start = '2001-01-01T00:00:00.000'; % info.timestr ignored.
%   opts.td.dt    = 1000;                      % 1-second cadence.

opts.td.window = struct();
    % Note: Same window applied to input and output
    opts.td.window.function = ''; 
    opts.td.window.functionstr = 'none';
    opts.td.window.functionargs = {};
    opts.td.window.width = NaN;   % Segment width.
    opts.td.window.shift = NaN;   % Shift amount to form new segment.
    opts.td.window.plot = [0,0,0];
    opts.td.window.loglevel = 0;

    % Example of applying a time domain window function to each segment.
    if 0
        % See 'help rectwin' for list of available MATLAB functions that
        % can be passed to tdwindow() (which is part of this package).
        opts.td.window.function = @tdwindow; 
        opts.td.window.functionstr = 'Rectangular';
        opts.td.window.functionargs = {@rectwin};
    end

opts.td.prewhiten = struct();
    % Note: Same prewhitening filter applied to input and output
    opts.td.prewhiten.function = '';
    opts.td.prewhiten.plot = [0,0,0];
    opts.td.prewhiten.loglevel = 0;
    
    % Example of prewhitening.
    if 0
        opts.td.prewhiten.function = @tdprewhiten;
        opts.td.prewhiten.functionstr = 'First difference';
        opts.td.prewhiten.functionargs = {'diff'};
    end

opts.fd.evalfreq = struct();
    opts.fd.evalfreq.plot = [0,0,0];
    opts.fd.evalfreq.loglevel = 0;

    % See evalfreq_demo.m for evalfreq() examples.
    opts.fd.evalfreq.function = @evalfreq;   
    opts.fd.evalfreq.functionargs = {7, 'logarithmic'};
    opts.fd.evalfreq.functionstr  = ...
                                    sprintf('%d frequencies per decade',...
                                    opts.fd.evalfreq.functionargs{1});

opts.fd.window = struct();
    opts.fd.window.function = @rectwin; 
    opts.fd.window.functionstr = 'rectangular';
    opts.fd.window.loglevel = 0; 
    opts.fd.window.plot = 0; 

opts.fd.stack = struct();
    opts.fd.stack.average = struct();
        opts.fd.stack.average.function = @transferfnAverage;
        opts.fd.stack.average.functionstr = '';
        opts.fd.stack.average.functionargs = {};

opts.fd.interpolation = struct();
    opts.fd.interpolation.function = @zinterp;
    opts.fd.interpolation.functionstr = 'zinterp()';
    opts.fd.interpolation.functionargs = {...
                    struct('loglevel',0,...
                           'interp1args',{{'linear',0}})};
    
opts.fd.regression = struct();
    opts.fd.regression.function = @ols_regress;
    opts.fd.regression.functionstr = 'OLS using regress() function';
    opts.fd.regression.functionargs = {struct('realvalued',0,'loglevel',0)};
    opts.fd.regression.plot = 0;
    opts.fd.regression.loglevel = 0;

    %opts.fd.regression.function = @ols_analytic;
    %opts.fd.regression.functionstr = 'OLS using analytic formula';
    %opts.fd.regression.functionargs = {};

    %opts.fd.regression.function = @robust_robustfit;
    %opts.fd.regression.functionstr = 'Robust regression using robustfit()';
    %opts.fd.regression.functionargs = {[],[],'off'};

    %opts.fd.regression.function = @robust_v1;
    %opts.fd.regression.functionstr = 'Robust regression using robust_v1()';
        %ropts = struct();
        %ropts.weightfn = 'huber';
        %ropts.stepmax = 50;
        %ropts.zeps = sqrt(eps);
        %ropts.hardcut = Inf; % 2.8
        %ropts.snstop = 1000;
        %ropts.verbose = 0;    
        %opts.fd.regression.functionargs = {ropts};

if os == 0
    % When no noise, should get exact TF used to generate the output data
    % (within limits of numerical precision).
    opts.description = 'OLS and $\Delta_f = 0$, $N_f = 0$.';
    opts.fd.evalfreq.functionargs = {[1,0], 'linear'};
    opts.fd.evalfreq.functionstr  = '1 DFT point per freq. band';
elseif os == 1 || nargs == 0
    opts.description = 'OLS and 7 pts/decade';
elseif os == 2
    opts.description = 'yulewalker(10) prewhiten';
    opts.td.prewhiten.method = 'yulewalker';
    opts.td.prewhiten.methodstr = 'yulewalker10';
    opts.td.prewhiten.options = 10;
elseif os == 3
    opts.description = 'Parzen window in FD';
    opts.fd.window.function = @parzenwin; 
    opts.fd.window.functionstr = 'parzen';
elseif os == 4
    opts.description = 'Robust regression';
    opts.fd.regression.method = 'robust_robustfit'; 
    opts.fd.regression.methodstr = 'Robust regression using robustfit()';
elseif os == 5
    opts.description = 'PCA rotation';
    opts.td.transform = 'pca';
    %opts.td.transform.methodstr = 'pca';
elseif os == 6
    opts.description = '1.5 day window';
    opts.td.window.width = 3600*48;
    opts.td.window.shift = 3600*24;
else
    error('Invalid option set number');
end
