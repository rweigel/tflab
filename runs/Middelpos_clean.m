function [B,E,t,infile,outfile] = Middelpos_clean(start,stop,rundir)

mat_raw = fullfile(rundir,'measurements-raw.mat');
mat_cleaned = fullfile(rundir,'measurements-cleaned.mat');

if ~exist(rundir,'dir')
    mkdir(rundir);
end
if exist(mat_cleaned,'file')
    fprintf('Reading: %s\n',mat_cleaned);
    load(mat_cleaned)
    return
end

if ~exist(mat_raw,'file')
    datadir = fullfile(rundir,'..','measurements');
    addpath(fullfile(scriptdir(),'readers'));
    if startsWith("2012", start) && startsWith("2012", stop)
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

E = despike(E,3,[-1,5]);

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
        E(I,:) = 0;
    end
end

%% Remove mean, interpolate over NaNs, pad E if NaNs at start or end
[B,E] = removemean(B,E);

infile = mat_raw;
outfile = mat_cleaned;
fprintf('Writing: %s\n',mat_cleaned);
save(mat_cleaned,'B','E','t','infile','outfile')
