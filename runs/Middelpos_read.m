%function [B,E,t] = Middelpos_read()

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

name = 'Dealesville';
ext = 't06';
inpath = [scriptpath,'data/Dealesville/measurements'];
outpath = [scriptpath,'data/Dealesville'];

if 1
name = 'Middelpos';
ext = 't82';
inpath = [scriptpath,'/data/Middelpos/measurements']; 
outpath = [scriptpath,'/data/Middelpos']; 
end

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
matfile = [outpath,'/',fname]; 
save(matfile,'B','E','t');
fprintf('Wrote: %s\n',matfile);
