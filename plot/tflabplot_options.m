function opts = tflabplot_options(S,opts,plotfun)

if iscell(S)
    nOut = size(S{1}.Out,2);
    nIn = size(S{1}.In,2);
    filestr = S{1}.Options.filestr;
    const_term = S{1}.Options.fd.regression.const_term;
else
    nOut = size(S.Out,2);
    nIn = size(S.In,2);
    filestr = S.Options.filestr;
    const_term = S.Options.fd.regression.const_term;
end

assert(isstruct(opts),'opts must be a struct.')

dopts = struct(); % Defalt opts

dopts.title = '';

% Print options passed to tflab's figsave()
dopts.print = 0;
dopts.printOptions.printName = plotfun;
dopts.printOptions.printDir = '';
dopts.printOptions.printFormats = {'pdf'};
dopts.printOptions.export_fig = {};
dopts.Position = [0 0 600 600];

% Line options passed to plot()
dopts.line = {'marker', '.', 'markersize', 15, 'linestyle', 'none'};

% Legend options passed to legend()
dopts.legend = {'Location', 'NorthWest', 'Orientation', 'Horizontal'};

% Options passed to subplot()
%                      [left bottom width height]
dopts.PositionTop    = [0.130 0.540 0.775 0.400];
dopts.PositionBottom = [0.130 0.110 0.775 0.400];

if strcmp(plotfun,'tsplot')
    if ~isfield(opts,'type') || isempty(opts.type)
        opts.type = 'original';
    end
    dopts.printOptions.printName = [dopts.printOptions.printName,'-',opts.type];
end

if strcmp(plotfun,'dftplot')
    if ~isfield(opts,'type') || isempty(opts.type)
        opts.type = 'original-raw-magnitudes';
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
    dopts.printOptions.printName = [dopts.printOptions.printName,'-',opts.type];
end

if strcmp(plotfun,'zplot')
    dopts.unwrap = 0;
    if ~isfield(opts,'type') || isempty(opts.type)
        opts.type = 1;
    end
    switch opts.type
        case 1
            dopts.printOptions.printName = [dopts.printOptions.printName,'-magnitude_phase'];
        case 2
            dopts.printOptions.printName = [dopts.printOptions.printName,'-rhoa_phase'];
        case 3
            dopts.printOptions.printName = [dopts.printOptions.printName,'-real_imaginary'];
    end
end

if strcmp(plotfun,'qqplot')
    if ~isfield(opts,'type') || isempty(opts.type)
        opts.type = 'standard';
    end
    % Determined by manually changing the figure size in the GUI and using
    % get(gcf,'Position').
    %dopts.Position = [0 0 391 625];
end

dopts.printOptions.printName = [dopts.printOptions.printName,'-',filestr];

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
        if const_term == 1
            dopts.zstrs = {'Z_{xx}','Z_{xy}','Z_{x0}', 'Z_{yx}','Z_{yy}', 'Z_{y0}'};
            dopts.rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{x0}','\rho^a_{yx}','\rho^a_{yy}','\rho^a_{y0}'};
            dopts.phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{x0}','\phi_{yx}','\phi_{yy}','\phi_{y0}'};
        end
    end
end

if any(strcmp(plotfun,{'dftplot','zplot','snplot'}))
    dopts.vs_period = 1;
    if iscell(S)
        for s = 1:length(S)
            freqsfs(s) = S{s}.Metadata.freqsf;
            T = 1./(S{s}.fe*freqsfs(s));
            Tmaxes(s) = max(T(isfinite(T)));
        end
        freqsf = min(freqsfs);
        dopts.period_range = [2/freqsf, max(Tmaxes)];
    else
        % TODO: Duplicate code.
        freqsf = S.Metadata.freqsf;
        T = 1./S.fe;
        Tmax = max(T(isfinite(T)));
        dopts.period_range = [2, Tmax]*(1/freqsf);
    end
    dopts.frequency_range = [0, 0.5]*freqsf;
end

if iscell(S)
    if isfield(S{1},'Metadata')
        dopts.instr = namestrs_(S{1}.Metadata.instr, nIn);
        dopts.outstr = namestrs_(S{1}.Metadata.outstr, nOut);
    end
else
    if isfield(S,'Metadata')
        dopts.instr = namestrs_(S.Metadata.instr, nIn);
        dopts.outstr = namestrs_(S.Metadata.outstr, nOut);
    end
end

opts = updateFields(opts,dopts);

end

% Replace default options with given options in opts
function opts = updateFields(opts, dopts)
    fns = fieldnames(dopts);
    for i = 1:length(fns)
        if ~isstruct(dopts.(fns{i}))
            if ~isfield(opts, fns{i})
                opts.(fns{i}) = dopts.(fns{i});
            end
        else
            if ~isfield(opts,fns{i})
                opts.(fns{i}) = dopts.(fns{i});
            else
                opts.(fns{i}) = updateFields(opts.(fns{i}),dopts.(fns{i}));
            end
        end
    end
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
