function h = tflab_title(tf, opts, plottype)

if ~isstruct(tf)
    return;
end

if ~isempty(opts.title)
    h = title(opts.title);
    return;
end

spacer = ' $|$ ';

tfo = tf.Options;
tfoi = tf.Options.info;

ts = '';
if ~isempty(tfoi.stationid) && isempty(tfoi.chainid)
    ts = sprintf('Site: %s',tfoi.stationid);            
end
if ~isempty(tfoi.stationid) && ~isempty(tfoi.chainid)
    ts = sprintf('%s/%s',tfoi.chainid, tfoi.stationid);
end
if ~isempty(ts)
    ts = [ts, spacer];
end    

if strcmp(plottype, 'ts')
    if strcmp(opts.type,'windowed') && ~isempty(tfo.td.window.function)
        ts = sprintf('%s%s windowed',...
                     ts,tfo.td.window.functionstr);
    end
    if strcmp(opts.type,'prewhitened') && ~isempty(tfo.td.prewhiten.function)
        ts = sprintf('%s%s prewhitened',...
                     ts,tfo.td.prewhiten.functionstr);
    end
    if strcmp(opts.type,'error')
        ts = sprintf('%s%s',ts,tfo.description);
    end
end

if strcmp(plottype, 'psd')
    if strcmp(opts.type,'zeropadded') && isfield(tf,'Zeropad')
        ts = sprintf('%sPadded with %d zeros',ts,tf.Options.td.zeropad);
    end    
    if strcmp(opts.type,'windowed') && isfield(tf,'Window')
        ts = [tf.Options.td.window.functionstr, '-windowed (in TD)'];
    end
    if strcmp(opts.type,'error')
        ts = sprintf('%s%s',ts,tfo.description);
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
