function [B,E,t,outfile,infiles] = LEMI_read(name, ext, inpath, outpath)

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

dlist = dir(inpath);

% Create .mat version of data in all text files
Data = [];
first = '';
for i = 1:length(dlist)
    fname = [dlist(i).folder,'/',dlist(i).name];
    fnamemat = [dlist(i).folder,'/',dlist(i).name,'.mat'];
    if ~endsWith(dlist(i).name,ext)
        continue;
    end
    try
        fprintf('Reading: %s\n',fname);
        data = load(fname);
    catch ex
        fprintf('Failed read. Skipping.\n');
        ex
    end
    infiles{i} = fname;
    fprintf('Writing: %s\n',fnamemat);
    save(fnamemat,'data');
    last = dlist(i).name(1:end-8);
    if isempty(first)
      first = dlist(i).name(1:end-8);
    end
    Data = [Data; data];
end

B = Data(:,7:9);
E = Data(:,12:13);
t = datenum(Data(:,1:6));

fname = sprintf('%s_%s-%s_raw.mat',name,first,last);
outfile = [outpath,'/',fname]; 
save(outfile,'B','E','t','outfile','infiles');
fprintf('Wrote: %s\n',outfile);
