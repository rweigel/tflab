function [B,E,t,infile,outfile] = Middelpos_clean(start,stop,rundir,usecache)

if nargin < 4
    usecache = 1;
end

mat_raw = fullfile(rundir,'measurements-raw.mat');
mat_cleaned = fullfile(rundir,'measurements-cleaned.mat');

logfile = fopen(fullfile(rundir,'Middelpos-cleaned.log'),'w');

if ~exist(rundir,'dir')
    mkdir(rundir);
end
if exist(mat_cleaned,'file') && usecache
    fprintf('Reading: %s\n',mat_cleaned);
    load(mat_cleaned);
    return
end

if ~exist(mat_raw,'file')
    datadir = fullfile(rundir,'..','measurements');
    if ~exist(datadir,'dir')
        error('Data directory not found: %s',datadir)
    end
    addpath(fullfile(scriptdir(),'readers'));
    if startsWith(start, "2012") && startsWith(stop, "2012")
        [B,E,t] = LEMI_read(datadir,'t82',start,stop,[50,50]);
    elseif str2num(start(1:4)) >= 2017
        [B,E,t] = LEMI_read(datadir,'t81',start,stop,[50,50]);
    else
        error('Time range not handled.')
    end
    fprintf('Writing: %s\n',mat_cleaned);
    save(mat_raw,'B','E','t');
else
    fprintf('Reading: %s\n',mat_raw);
    load(mat_raw);
end

msg = 'Despiking E\n';
logmsg(msg);
fprintf(logfile,msg);
E = despike(E,0.1,[1,5],logfile);

% TODO: Report on largest gap
ti = [1:size(B,1)]';
for i = 1:size(B,2)
    tg = find(~isnan(B(:,i)));
    B(:,i) = interp1(tg,B(tg,i),ti);
end
for i = 1:size(E,2)
    tg = find(~isnan(E(:,i)));
    E(:,i) = interp1(tg,E(tg,i),ti);
end

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
