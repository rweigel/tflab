function [B,E,t] = Middelpos_data()

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Dir of .t82 files
datapath = [scriptpath,'/data/Middelpos']; 

% Content of t82 files
matfile = [datapath,'/Middelpos-measurements.mat']; 

if ~exist(matfile,'file')
    dname = [datapath,'/measurements'];
    do = datenum(2012,7,12);
    df = datenum(2012,9,4);
    D = [];
    for d = do:df
        [y,m,d] = datevec(d);
        fname = sprintf('%s/MT_Middelpos_%d%02d%02d.t82',dname,y,m,d);
        fprintf('Reading %s\n',fname);
        tmp = load(fname);
        D = [D;tmp];
    end
    B = D(:,7:9);
    E = D(:,12:13);
    t = datenum(D(:,1:6));
    save(matfile,'D','B','E','t');
else
    load(matfile);
end
