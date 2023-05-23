function [B,E,t,infile,outfile] = Middelpos_clean()

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

mat_raw = [scriptpath,'/data/Middelpos/Middelpos_20120712-20121107_raw.mat']; 
mat_clean = [scriptpath,'/data/Middelpos/Middelpos_20120712-20121107_cleaned.mat']; 

if exist(mat_clean,'file')
    fprintf('Reading: %s\n',mat_clean);
    load(mat_clean)
    return
end

if ~exist(mat_raw,'file')
    LEMI_read(name, ext, inpath, outpath);
end

fprintf('Reading: %s\n',mat_raw);
load(mat_raw)

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
        logmsg(sprintf('Set %d leading or trailing NaNs in column %d of E to zero\n',length(I),i))
        E(I,:) = 0;
    end
end

%% Remove mean, interpolate over NaNs, pad E if NaNs at start or end
[B,E] = removemean(B,E);

infile = mat_raw;
outfile = mat_clean;
fprintf('Writing: %s\n',mat_clean);
save(mat_clean,'E','B','t','infile','outfile')
