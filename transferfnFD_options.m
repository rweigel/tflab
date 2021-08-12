function opts = transferfnFD_options(os,iopts)
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
    opts.info.instr = 'In'; % Cell array or string. 
    % Or {'$B_x$', ...} (1 cell element per column in In)
    
    opts.info.outstr = 'Out'; % Cell array or string. 
    % Or {'$E_x$', ...} (1 cell element per column in Out)
    
    opts.info.timestart = ''; 
    opts.info.timedelta = 1; % Measurement cadence

    % Timestart is a time string of the form
    % 'yyyy-mm-ddTHH:MM:SS.FFF'.
    % Example: 
    %   opts.td.timestart = '2001-01-01T00:00:00.000';
    
    %opts.info.inunit= 'nT';
    %opts.info.outunit= 'mV/km';
    %opts.info.timeunit = 's';

    opts.info.inunit= '';
    opts.info.outunit= '';
    opts.info.timeunit = '';
    opts.info.stationid = '';
    opts.info.chainid = '';
    
opts.transferfnFD = struct();
    opts.transferfnFD.loglevel = 1;
    % Elements of arrays for plot.* are 0 or 1 and
    % [showplot, savepng, savepdf] 
    opts.transferfnFD.plot = struct();
        opts.transferfnFD.plot.timeseries = [0,0,0];
        opts.transferfnFD.plot.spectrum = [0,0,0];
        opts.transferfnFD.plot.Z = [0,0,0];
        opts.transferfnFD.plot.H = [0,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time domain options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% # of points at start and end to trim before computing metrics
% (pe/cc/mse/sn/coherence)
opts.td.Ntrim = NaN;

opts.td.detrend.function = struct();
    %opts.td.detrend.function = @removemean;
    opts.td.detrend.function = '';
    opts.td.detrend.functionstr = '';  % Optional descriptive name
    opts.td.detrend.functionargs = {}; % Arguments after first argument to fn.

% Dimensionless start; ignored if time array passed to transferfnFD.
opts.td.start = 1; 

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
        % Default is equivalent to
        opts.td.window.function = @tdwindow; 
        opts.td.window.functionstr = 'Rectangular';
        opts.td.window.functionargs = {@rectwin};
        %opts.td.window.functionstr = 'Parzen';        
        %opts.td.window.functionargs = {@parzenwin};        
    end

opts.td.prewhiten = struct();
    % Note: Same prewhitening filter applied to input and output
    opts.td.prewhiten.function = '';
    opts.td.prewhiten.plot = [0,0,0];
    opts.td.prewhiten.loglevel = 0;
    
    % Example of prewhitening.
    if 0
        opts.td.prewhiten.function = @prewhiten;
        opts.td.prewhiten.functionstr = 'First difference';
        opts.td.prewhiten.functionargs = {'diff'};
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.fd.progam = struct();
    opts.fd.program.name = 'transferfnFD';
    %opts.fd.progam.name = 'lemimt'; % Use lemimt program. All fd options below are ignored.
    opts.fd.program.options = '';

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

% When stack.average.function = @transferfnAverage is used, ffts in each
% time window are computed and the ffts in a frequency band are used for a
% regression to find Z for that time window and frequency band. Then, the
% Zs for each time window are averaged.
% When opts.fd.stack.average.function == '', ffts in each time window are
% computed and the ffts in a frequency band for all windows are combined
% and used to perform one regression to compute one Z in that frequency
% band.
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
    opts.fd.regression.functionargs = {struct('realvalued',0,'loglevel',1)};
    opts.fd.regression.plot = 0;
    opts.fd.regression.loglevel = 1;

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
    opts.description = 'OLS and $\Delta {f_e} = 1$, $\Delta w = 0$.';
    opts.fd.evalfreq.functionargs = {[1,0], 'linear'};
    opts.fd.evalfreq.functionstr  = '1 DFT point per freq. band';
elseif os == 1 || nargs == 0
    opts.description = 'OLS and $\sim$7 pts/decade';
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

if nargin > 1
    opts = options(opts,iopts);
end

function opts = options(opts,iopts)
    % Overwrite default options if options given
    fns = fieldnames(iopts);
    for i = 1:length(fns)
        if isfield(opts,fns{i})
            if isstruct(iopts.(fns{i}))
                opts.(fns{i}) = options(opts.(fns{i}),iopts.(fns{i}));
            else
                opts.(fns{i}) = iopts.(fns{i});
            end
        else
            warning(sprintf('Option %s is not a valid option',fns{i}));
        end
    end
end

end