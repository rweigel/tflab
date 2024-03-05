function [B,E,t,infiles] = LEMI_read(inpath, ext, start, stop)
% LEMI_read - Read LEMI data from text files
%
%   [B,E,t,infiles] = LEMI_read(inpath, ext, start, stop)
%
%   Assumes inpath contains LEMI text files with names that start with
%   YYYYMMMDD and end with .ext.
%
%   start and stop are strings that specify the day range of files to read
%   (inclusive), and must have the format YYYYMMMDD.
%
%   Creates .mat version of data in all text files in the same directory.

dlist = dir(inpath);

Data = [];
first = '';
for i = 1:length(dlist)
    fname = [dlist(i).folder,'/',dlist(i).name];
    fnamemat = [dlist(i).folder,'/',dlist(i).name,'.mat'];
    if ~endsWith(dlist(i).name,ext)
        continue
    end
    if dlist(i).name(1:8) < start || (dlist(i).name(1:8) > stop)
        continue
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
