function [B,E,t,infile,outfile] = Middelpos_clean(start,stop,basedir)

mat_raw = sprintf('Middelpos_%s-%s_raw.mat',start,stop);
mat_clean = sprintf('Middelpos_%s-%s_cleaned.mat',start,stop);

mat_raw = fullfile(basedir,'data','Middelpos',mat_raw);
mat_clean = fullfile(basedir,'data','Middelpos',mat_clean);

if exist(mat_clean,'file')
    fprintf('Reading: %s\n',mat_clean);
    load(mat_clean)
    return
end

if ~exist(mat_raw,'file')
    [B,E,t] = LEMI_read(inpath, ext, start, stop);
    fprintf('Writing: %s\n',mat_clean);
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
outfile = mat_clean;
fprintf('Writing: %s\n',mat_clean);
save(mat_clean,'B','E','t','infile','outfile')
