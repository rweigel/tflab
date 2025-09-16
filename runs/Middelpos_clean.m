function [B,E,t,infile,outfile] = Middelpos_clean(start,stop,rundir,usecache)

if nargin < 4
    usecache = 1;
end

mat_raw = fullfile(rundir,'measurements-raw.mat');
mat_cleaned = fullfile(rundir,'measurements-cleaned.mat');

if ~exist(rundir,'dir')
    mkdir(rundir);
end

logfile = fopen(fullfile(rundir,'Middelpos-cleaned.log'),'w');

if exist(mat_cleaned,'file') && usecache
    fprintf('Reading: %s\n',mat_cleaned);
    load(mat_cleaned);
    return
end

if exist(mat_raw,'file') && usecache
    fprintf('Reading: %s\n',mat_raw);
    load(mat_raw);
else
    datadir = fullfile(rundir,'..','measurements');
    if ~exist(datadir,'dir')
        error('Data directory not found: %s',datadir)
    end
    addpath(fullfile(scriptdir(),'readers'));
    % LEMI instrument was configured with lengths = 1 m. Here we correct
    % to set distance between ends of electric field probes to be 100 m.
    if startsWith(start, "2012") && startsWith(stop, "2012")
        [B,E,t] = LEMI_read(datadir,'t82',start,stop,[100,100]);
    elseif str2num(start(1:4)) >= 2017
        [B,E,t] = LEMI_read(datadir,'t81',start,stop,[100,100]);
    else
        error('Time range not handled.')
    end
    size(B)
    fprintf('Writing: %s\n',mat_cleaned);
    save(mat_raw,'B','E','t');
end

msg = 'Despiking E\n';
logmsg(msg);
fprintf(logfile,msg);
E = despike(E,0.1,[1,5],logfile);

% TODO: Report on largest gap
ts  = round(86400*(t-t(1))); % Time in seconds since start
ti = ts(1):ts(end); % Interpolation grid
msg = sprintf('Max dt = %d [s]\\n', max(diff(ts)));
logmsg(msg)
fprintf(logfile, msg);

for i = 1:size(B,2)
    tg = find(~isnan(B(:,i)));
    msg = sprintf('B(:,%d) has %d NaNs\\n', i, length(t)-length(tg));
    logmsg(msg)
    fprintf(logfile, msg);
    B(:,i) = interp1(ts(tg),B(tg,i),ti);
end
for i = 1:size(E,2)
    tg = find(~isnan(E(:,i)));
    msg = sprintf('E(:,%d) has %d NaNs\\n', i, length(t)-length(tg));
    logmsg(msg)
    fprintf(logfile, msg);
    E(:,i) = interp1(ts(tg),E(tg,i),ti);
end
t = t(1) + ti/86400; % Interpolation grid time in datenum

for i = 1:size(E,2)
    I = find(isnan(E(:,i)));
    if ~isempty(I)
        msg = 'Set %d leading or trailing NaNs in column %d of E to zero\n';
        logmsg(msg,length(I),i);
        fprintf(logfile, '%s', msg);
        E(I,:) = 0;
    end
end

[B,E] = removemean(B,E);

infile = mat_raw;
outfile = mat_cleaned;
logmsg('Writing: %s\n',mat_cleaned);
save(mat_cleaned,'B','E','t','infile','outfile')
fclose(logfile);
