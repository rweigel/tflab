function h = tflab_title(tf, opts, plottype)

if ~isempty(opts.title)
    h = title(opts.title);
    return;
end

tfo = tf.Options;
tfoi = tf.Options.info;
ts = '';

if strcmp(plottype, 'ts')
    if ~isempty(tfoi.stationid) && isempty(tfoi.chainid)
        ts = sprintf('Site: %s',tfoi.stationid);            
    end
    if ~isempty(tfoi.stationid) && ~isempty(tfoi.chainid)
        ts = sprintf('%s/%s',tfoi.chainid, tfoi.stationid);
    end
    if strcmp(opts.type,'windowed') && ~isempty(tfo.td.window.function)
        ts = sprintf('%s | %s-windowed',...
                     ts,tfo.td.window.functionstr);
    end
    if strcmp(opts.type,'prewhitened') && ~isempty(tfo.td.prewhiten.function)
        ts = sprintf('%s | %s prewhitened',...
                     ts,tfo.td.prewhiten.functionstr);
    end
end

if strcmp(plottype, 'psd')
    if isstruct(tf)
        if strcmp(opts.type,'error')
            ts = tfo.description;
        end
    end
end

if any(strcmp(plottype, {'z','sn'})) && isstruct(tf)
    ts = tfo.description;
    if ~isempty(tfoi.stationid) && isempty(tfoi.chainid)
        ts = sprintf('%s | %s',tfoi.stationid, ts);            
    end
    if ~isempty(tfoi.stationid) && ~isempty(tfoi.chainid)
        ts = sprintf('%s/%s | %s',tfoi.chainid, tfoi.stationid, ts);
    end
end

h = [];
if ~isempty(ts)
    h = title(ts);
end
