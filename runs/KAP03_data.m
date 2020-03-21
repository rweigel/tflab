function [B,E] = KAP_data(id)

addpath([fileparts(mfilename('fullpath')),'/../misc']); % logmsg.m

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Dir of data file
datapath = [scriptpath,'/data/KAP03']; 

% Content of data file
matfile = [datapath,'/',id,'-measurements.mat']; 

if ~exist(matfile,'file')
    fname = [datapath,'/measurements/',lower(id),'as.ts']
    logmsg(sprintf('Reading %s\n',fname));
    fid = fopen(fname);
    data = textscan(fid,'%f %f %f %f %f','CollectOutput',1,'HeaderLines',113) ;
    fclose(fid);

    % Replace missing data value with NaN
    data{1}(data{1} == 1.00000003e32) = nan;

    B = data{1}(:,1:3);
    E = data{1}(:,4:5);

    B(B == 1.00000003e32) = nan;
    E(abs(E)>30) = nan; % Obvious spikes

    logmsg('Saving %s\n',matfile);
    save(matfile,'E','B');
else
    logmsg(sprintf('Reading %s\n',matfile));
    load(matfile);
end
