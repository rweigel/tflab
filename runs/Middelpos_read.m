%function [B,E,t] = Middelpos_read()

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Dir of .t82 files
datapath = [scriptpath,'/data/Middelpos']; 

datapath = [scriptpath,'/data/Middelpos/measurements']; 
dlist = dir(datapath);

% Create .mat version of data in .t82 files
for i = 1:length(dlist)
    fname = [dlist(i).folder,'/',dlist(i).name];
    if endsWith(dlist(i).name,'t82') && ~exist(fname,'file')
        fprintf('Reading: %s\n',fname);
        data = load(fname);
        save([fname,'.mat'],'data');
        fprintf('Wrote:   %s\n',[fname,'.mat']);                
    end
end

% Read and concatenate .mat files
Data = [];
first = '';
for i = 1:length(dlist)
    if endsWith(dlist(i).name,'mat')
        last = dlist(i).name(1:end-8);
        if isempty(first)
          first = dlist(i).name(1:end-8);
        end
        fname = [dlist(i).folder,'/',dlist(i).name];
        fprintf('Reading: %s\n',fname);                
        data = load(fname);
        Data = [Data; data.data];
    end
end

B = Data(:,7:9);
E = Data(:,12:13);
t = datenum(Data(:,1:6));

datapath = [scriptpath,'/data/Middelpos/']; 
fname = sprintf('Middelpos_%s-%s_raw.mat',first,last);
matfile = [datapath,'/',fname]; 
save(matfile,'B','E','t');
fprintf('Wrote: %s\n',matfile);
