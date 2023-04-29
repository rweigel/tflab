function opts = tflab_options(os,iopts)
%TFLAB_OPTIONS - Return options for tflab().
%
%  opts = TFLAB_OPTIONS() returns default options.
%
%  opts = TFLAB_OPTIONS(os), where os is an option set.
%
%  opts = TFLAB_OPTIONS(os, opts1) returns opts with values in opts1 by
%  calling updatedefaults. opts1 does not need to contain all options
%  returned by TFLAB_OPTIONS.
% 
%  See also TFLAB, TFLAB_OPTIONS, UPDATEDEFAULTS.

if nargin == 0
    os = 0;
end

opts = struct();

opts.filestr = sprintf('tflab_options-%d',os);

opts.tflab = struct();
    opts.tflab.loglevel = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot label options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.info = struct();

    opts.info.instr = 'In'; % Cell array or string. 
    % Or ,{'$B_x$', ...} (1 cell element per column in In)
    
    opts.info.outstr = 'Out'; % Cell array or string. 
    % Or, {'$E_x$', ...} (1 cell element per column in Out)
    
    % Timestart is a time string of the form
    %   'yyyy-mm-ddTHH:MM:SS.FFF'.
    % Example: 
    %   opts.td.timestart = '2001-01-01T00:00:00.000';
    opts.info.timestart = ''; 
    
    % Measurement cadence. Specify units using opt.info.timeunit.  
    opts.info.timedelta = 1; 
    
    %opts.info.inunit= 'nT';
    %opts.info.outunit= 'mV/km';
    %opts.info.timeunit = 's';

    opts.info.inunit    = '';
    opts.info.outunit   = '';
    opts.info.timeunit  = '';   % Unit for timedelta
    opts.info.chainid   = '';   % Usualy over-arching project name
    opts.info.stationid = '';   % Usually an abbreviation, e.g., VAQ58

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time domain options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% # of points at start and end to trim before computing metrics
% (pe/cc/mse/sn/coherence). NaN => no trim.
opts.td.Ntrim = 60*30;

% Dimensionless start time; ignored if time array passed to tflab.
opts.td.start = 1;

% Number of zeros added to the end of all time series prior to computing
% DFT. NaN => no pad.
opts.td.zeropad = NaN;

opts.td.detrend.function = struct();
    %opts.td.detrend.function = @removemean;
    opts.td.detrend.function = '';
    opts.td.detrend.functionstr = '';  % Optional descriptive name
    opts.td.detrend.functionargs = {}; % Arguments after first argument to fn.

opts.td.window = struct();
    % Note: Same window applied to input and output
    opts.td.window.function = ''; 
    opts.td.window.functionstr = 'none';
    opts.td.window.functionargs = {};
    opts.td.window.width = NaN;   % Segment width.
    opts.td.window.shift = NaN;   % Shift amount to form new segment.
    opts.td.window.loglevel = 0;

    % Example of applying a time domain window function to each segment.
    if 0
        % See https://www.mathworks.com/help/signal/ug/windows.html
        % for list of available MATLAB functions that can be passed to
        % tdwindow() (which is function in this package).
        % Default is equivalent to
        %opts.td.window.function = @tdwindow; 
        %opts.td.window.functionstr = 'Rectangular';
        %opts.td.window.functionargs = {@rectwin};
    end

opts.td.prewhiten = struct();
    % Note: Same prewhitening filter applied to input and output
    opts.td.prewhiten.function = '';
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

opts.fd.program = struct();
    opts.fd.program.name = 'tflab';
    % Use lemimt program. All fd options below are ignored.
    %opts.fd.program.name = 'lemimt';
    opts.fd.program.options = '';

opts.fd.evalfreq = struct();
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
    opts.fd.regression.function = @regress_ols;
    opts.fd.regression.functionstr = 'OLS using backslash function';
    opts.fd.regression.functionargs = {'backslash'};
    %opts.fd.regression.function = @ols_regress;    
    %opts.fd.regression.functionstr = 'OLS using regress() function';
    %opts.fd.regression.functionargs = {struct('realvalued',0,'loglevel',0)};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.plot = struct();
    opts.plot.figprep = @figprep; % Function to call prior to plotting
    opts.plot.figsave = @figsave; % Function to call for printing

    opts.plot.tsplot.type = 'raw';
    opts.plot.tsplot.title = '';
    opts.plot.tsplot.print = 0;
    opts.plot.tsplot.printname = 'tsplot';
    opts.plot.tsplot.printdir = '';
    opts.plot.tsplot.printfmt = {'pdf'};
    % printfmt is one or more extensions allowed by export_fig()
    % e.g. {'pdf'} or {'svg','png'}.

if os == 0
    % When no noise, should get exact TF used to generate the output data
    % (within limits of numerical precision).
    fstr = '1 point per freq. band.';
    opts.description = ['OLS with ', fstr];
    opts.fd.evalfreq.functionargs = {[1,0], 'linear'};
    opts.fd.evalfreq.functionstr  = fstr;
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
    opts = updatedefaults(opts,iopts);
end