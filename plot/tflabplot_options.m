function opts = tflabplot_options(S,opts,dtype,plotfun)

dopts = struct(); % Defalt opts

dopts.title = '';

dopts.type = dtype;

prefix = struct('tsplot','ts',...
                'dftplot','dft',...
                'zplot','tf',...
                'snplot','sn',...
                'qqplot','qq');

% Print options passed to tflab's figsave()
dopts.print = 0;
dopts.printname = prefix.(plotfun);
dopts.printdir = '';
dopts.printfmt = {'pdf'};

% Line options passed to plot()
dopts.line = {'marker', '.', 'markersize', 15, 'linestyle', 'none'};

% Legend options passed to legend()
dopts.legend = {'Location', 'NorthWest', 'Orientation', 'Horizontal'};

% Options passed to subplot()
dopts.PositionTop = [0.1300 0.5400 0.7750 0.4];
dopts.PositionBottom = [0.1300 0.1100 0.7750 0.4];

if iscell(S)
    nOut = size(S{1}.Out,2);
    nIn = size(S{1}.In,2);
else
    nOut = size(S.Out,2);
    nIn = size(S.In,2);
end

if strcmp(plotfun,'dftplot')
    if isempty(opts.type)
        opts.type = 'original-raw-magphase';
    end
    tparts = split(opts.type,'-');
    if length(tparts) == 1
        tparts{2} = 'raw';
    end
    if length(tparts) == 2
        if strcmp(tparts{1},'error')
            tparts{3} = 'magphase';
        else
            tparts{3} = 'magnitudes';
        end
    end
    opts.type = join(tparts,'-');
    opts.type = opts.type{:};
end

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

if any(strcmp(plotfun,{'zplot','qqplot'}))
    if nOut == 1
        if nIn == 1
            dopts.zstrs = {'Z'};
            dopts.rhostrs = {'\rho^a'};
            dopts.phistrs = {'\phi'};
        end
        if nIn == 2
            dopts.zstrs = {'Z_{x}','Z_{y}'};
            dopts.rhostrs = {'\rho^a_{x}','\rho^a_{y}'};
            dopts.phistrs = {'\phi_{x}','\phi_{y}'};
        end
        if nIn == 3
            dopts.zstrs = {'Z_{x}','Z_{y}','Z_{z}'};
            dopts.rhostrs = {'\rho^a_{x}','\rho^a_{y}','\rho^a_{z}'};
            dopts.phistrs = {'\phi_{x}','\phi_{y}','\phi_{z}'};
        end
    end
    if nOut == 2 && nIn == 2
        dopts.zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
        dopts.rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{yx}','\rho^a_{yy}'};
        dopts.phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
    end    
end

if any(strcmp(plotfun,{'dftplot','zplot','snplot'}))
    dopts.vs_period = 1;
    if iscell(S)
        for s = 1:length(S)
            timedeltas(s) = S{s}.Metadata.timedelta;
            T = timedeltas(s)./S{s}.fe;
            Tmaxes(s) = max(T(~isnan(T)));
        end
        timedelta = min(timedeltas);
        dopts.period_range = [2*timedelta, max(Tmaxes)];
    else
        timedelta = S.Metadata.timedelta;
        T = 1./S.fe;
        Tmax = max(T(~isnan(T)));
        dopts.period_range = [2, Tmax]*timedelta;
    end
    dopts.frequency_range = [0, 0.5]/timedelta;     
end

if iscell(S)
    dopts.instr = namestrs_(S{1}.Metadata.instr, nIn);
    dopts.outstr = namestrs_(S{1}.Metadata.outstr, nOut);
else
    dopts.instr = namestrs_(S.Metadata.instr, nIn);
    dopts.outstr = namestrs_(S.Metadata.outstr, nOut);
end

% Replace default options with given options in opts
fns = fieldnames(dopts);
for i = 1:length(fns)
    if isfield(opts, fns{i})
        dopts.(fns{i}) = opts.(fns{i});
    end
end

opts = dopts;

end

function namestrs = namestrs_(namestrs, nc)
    if nc > 1
        if ~iscell(namestrs)
            for j = 1:nc
                namestrs_new{j} = sprintf('%s(:,%d)',namestrs,j);
            end
            namestrs = namestrs_new;
        else
            if length(namestrs) ~= nc
                error('Metadata.instr or Metadata.outstr does not have the correct number of elements.');
            end
        end
    else
        namestrs = {namestrs};
    end    
end
