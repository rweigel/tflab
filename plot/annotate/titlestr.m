function h = titlestr(tf, opts, plottype)

if ~isstruct(tf)
    return;
end

if ~isempty(opts.title) || ~isfield(tf,'Options')
    h = title(opts.title);
    return;
end

spacer = ' $|$ ';

tfom = tf.Metadata;
ts = '';
if ~isempty(tfom.stationid) && isempty(tfom.chainid)
    ts = sprintf('Site: %s',tfom.stationid);            
end
if ~isempty(tfom.stationid) && ~isempty(tfom.chainid)
    ts = sprintf('%s/%s',tfom.chainid, tfom.stationid);
end
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
        ftype = opts.type(1:end-2); % Remove "ed"
        if isfield(tfo,'td') && isfield(tfo.td, ftype) && ~isempty(tfo.td.(ftype).function)
            fdesc = tfo.td.(ftype).functionstr;
            if isempty(fdesc)
                ts = sprintf('%s%sed',ts,ftype);
            else
                ts = sprintf('%s%s %sed',ts,fdesc,ftype);
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

if any(strcmp(plottype, {'z','sn'})) && isstruct(tf)
    ts = sprintf('%s%s',ts,tfo.description);
end

if endsWith(ts,spacer)
    ts = ts(1:end-length(spacer));
end

h = [];
if ~isempty(ts)
    h = title(ts);
end
