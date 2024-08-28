function opts = tflabplot_options(S,opts,plotfun)

if ~iscell(S)
    S = {S};
end

filestr = S{1}.Options.filestr;
for i = 2:length(S)
    filestr = [filestr,';',S{i}.Options.filestr];
end

assert(isstruct(opts),'opts must be a struct.')

dopts = struct(); % Defalt opts

%dopts.title = '';

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

if isfield(S{1},'Metadata')
    for s = 1:length(S)
        nOut(s) = size(S{s}.Out,2);
        nIn(s) = size(S{s}.In,2);
        dopts.instr = namestrs_(S{1}.Metadata.instr, nIn);
        dopts.outstr = namestrs_(S{1}.Metadata.outstr, nOut);
    end
end

if any(strcmp(plotfun,{'zplot','qqplot'}))
    for s = 1:length(S)
        nOut(s) = size(S{s}.Out,2);
        nIn(s) = size(S{s}.In,2);
        const_term = S{s}.Options.fd.regression.const_term;
        if nOut(s) == 1
            if nIn(s) == 1
                dopts.zstrs{s} = {'Z'};
                dopts.rhostrs{s} = {'\rho^a'};
                dopts.phistrs{s} = {'\phi'};
            end
            if nIn(s) == 2
                dopts.zstrs{s} = {'Z_{x}','Z_{y}'};
                dopts.rhostrs{s} = {'\rho^a_{x}','\rho^a_{y}'};
                dopts.phistrs{s} = {'\phi_{x}','\phi_{y}'};
            end
            if nIn(s) == 3
                dopts.zstrs{s} = {'Z_{x}','Z_{y}','Z_{z}'};
                dopts.rhostrs{s} = {'\rho^a_{x}','\rho^a_{y}','\rho^a_{z}'};
                dopts.phistrs{s} = {'\phi_{x}','\phi_{y}','\phi_{z}'};
            end
            if const_term == 1
                dopts.zstrs{s} = {dopts.zstrs{:}, '\delta E'};
                dopts.rhostrs{s} = {dopts.rhostrs{:}, '\delta \rho'};
                dopts.phistrs{s} = {dopts.phistrs{:},'\delta \phi'};
            end
        end
        if nOut(s) == 2 && nIn(s) == 2
            for i = 1:2
                for j = 1:2
                    dopts.zstrs{s}{i,j} = sprintf('Z_{%d%d}',i,j);
                    dopts.rhostrs{s}{i,j} = sprintf('\\rho^a_{%d%d}',i,j);
                    dopts.phistrs{s}{i,j} = sprintf('\\phi_{%d%d}',i,j);
                end
                if const_term == 1
                    dopts.zstrs{s}{i,j+1} = sprintf('\\delta E_%d', i);
                    dopts.rhostrs{s}{i,j+1} = sprintf('\\delta \\rho^a_%d', i);
                    dopts.phistrs{s}{i,j+1} = sprintf('\\delta \\phi_%d', i);
                end
            end
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
        if ~iscell(namestrs)
            namestrs = {namestrs};
        end
    end
end
