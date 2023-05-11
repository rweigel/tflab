function opts = tflabplot_options(S,opts,dtype,plotfun)

% TODO: Can get plotfun from calling function.

dopts = struct();

dopts.title = '';

dopts.type = dtype;

prefix = struct('tsplot','ts','psdplot','psd','zplot','tf','snplot','sn');

% Print options passed to tflab's figsave()
dopts.print = 0;
dopts.printname = prefix.(plotfun);
dopts.printdir = '';
dopts.printfmt = {'pdf'};

% Line options passed to plot()
dopts.line = {'marker', '.', 'markersize', 20, 'linestyle', 'none'};

% Legend options passed to legend()
dopts.legend = {'Location', 'NorthWest', 'Orientation', 'Horizontal'};

% Options passed to subplot()
dopts.PositionTop = [0.1300 0.5400 0.7750 0.4];
dopts.PositionBottom = [0.1300 0.1100 0.7750 0.4];

if strcmp(plotfun,'zplot')
    dopts.unwrap = 0;
    switch dtype
        case 1
            dopts.printname = [dopts.printname,'_magnitude_phase'];
        case 2
            dopts.printname = [dopts.printname,'_rhoa_phase'];            
        case 3
            dopts.printname = [dopts.printname,'_real_imaginary'];            
    end
end

if any(strcmp(plotfun,{'psdplot','zplot','snplot'}))
    dopts.vs_period = 1;
    if iscell(S)
        for s = 1:length(S)
            lens(s) = size(S{s}.In,1);
        end
        dopts.period_range = [1, 2*max(lens)];
    else        
        dopts.period_range = [1, 2*size(S.In,1)];
    end
    dopts.frequency_range = [0, 0.5];
end

% Replace default options with given options in opts
fns = fieldnames(dopts);
for i = 1:length(fns)
    if isfield(opts, fns{i})
        dopts.(fns{i}) = opts.(fns{i});
    end
end

opts = dopts;
