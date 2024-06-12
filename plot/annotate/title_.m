function h = title_(tf, opts, plottype)

if ~isfield(opts,'title')
    return;
end

if length(opts.title) == 0
    return
end

if ~isempty(opts.title)
    h = title(opts.title);
    return;
end

if iscell(tf)
    for i = 1:length(tf)
        tsc{i} = titlestr_site(tf{i}.Metadata);
    end
    if length(unique(tsc)) == 1
        % Plot contains data from multiple tfs. If all titles are the same, 
        % use that title.
        h = title(tsc{1});
        return
    end
end

if ~isfield(tf,'Options')
    h = title(opts.title);
    return;
end

spacer = ' $|$ ';

ts = titlestr_site(tf.Metadata);
if ~isempty(ts)
    ts = [ts, spacer];
end

tfo = tf.Options;
if strcmp(plottype, 'ts')
    if strcmp(opts.type,'error')
        ts = sprintf('%s%s',ts,tfo.description);
    end
    if strcmp(opts.type,'zeropadded')
        ts = sprintf('%sPadded with %d zeros',ts,tf.Options.td.zeropad);
    else
        if strcmp(opts.type,'final')
            ts = sprintf('%sAfter preprocessing',ts);
        else
            ftype = opts.type(1:end-2); % Remove "ed"
            if isfield(tfo,'td') && isfield(tfo.td, ftype) && ~isempty(tfo.td.(ftype).function)
                fdesc = replace(tfo.td.(ftype).functionstr,'_','\_');
                if isempty(fdesc)
                    ts = sprintf('%s%sed',ts,ftype);
                else
                    ts = sprintf('%s$\\texttt{%s()}$ %sed',ts,fdesc,ftype);
                end
            end
        end
    end
end

if strcmp(plottype, 'dft')
    tparts = split(opts.type,'-');
    ftype = tparts{1}(1:end-2); % Remove "ed"
    if strcmp(tparts{1},'error')
        ts = sprintf('%s%s',ts,tfo.description);
    end
    if strcmp(tparts{1},'zeropadded')
        ts = sprintf('%sPadded with %d zeros',ts,tf.Options.td.zeropad);
    else
        if isfield(tfo, 'td') && isfield(tfo.td, ftype) && ~isempty(tfo.td.(ftype).function)
            fdesc = tf.Options.td.(ftype).functionstr;
            if isempty(fdesc)
                ts = sprintf('%s%sed',ts,ftype);
            else
                ts = sprintf('%s%s %sed',ts,fdesc,ftype);
            end
        end
    end
end

if any(strcmp(plottype, {'z','sn','qq'})) && isstruct(tf)
    ts = sprintf('%s%s',ts,tfo.description);
end

if endsWith(ts,spacer)
    ts = ts(1:end-length(spacer));
end

h = [];
if ~isempty(ts)
    h = title(ts);
end

end

function ts = titlestr_site(Metadata)
    ts = '';
    if ~isempty(Metadata.stationid) && isempty(Metadata.chainid)
        ts = sprintf('Site: %s',Metadata.stationid);
    end
    if ~isempty(Metadata.stationid) && ~isempty(Metadata.chainid)
        ts = sprintf('%s/%s',Metadata.chainid, Metadata.stationid);
    end
end