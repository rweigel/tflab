function [B,E,Header] = KAP_data(id)

addpath([fileparts(mfilename('fullpath')),'/../misc']); % logmsg.m

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Dir of data file
datapath = [scriptpath,'/data/KAP03']; 

% Content of data file
matfile = [datapath,'/',id,'-measurements.mat']; 

if 1 || ~exist(matfile,'file')
    fname = [datapath,'/measurements/',lower(id),'as.ts'];
    fid = fopen(fname);
    i = 1;
    while 1
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        headerlines{i} = line;
        if strcmp(line(1),'>')
            tmp = strsplit(line,':');
            if length(tmp) > 2
                % Line had multiple colons. Recover them.
                Header.(strtrim(tmp{1}(2:end))) = strtrim(strjoin(tmp(2:end),':'));
            else
                Header.(strtrim(tmp{1}(2:end))) = strtrim(tmp{2});
            end

        end
        if strncmp(line,'>INFO_END',9)
            break
        end
        i = i + 1;
    end
    fclose(fid);
    Header.headerlines = headerlines;

    logmsg(sprintf('Reading %s\n',fname));
    fid = fopen(fname);
    data = textscan(fid,'%f %f %f %f %f','CollectOutput',1,'HeaderLines',113) ;
    fclose(fid);

    MIS_DATA = str2num(Header.MIS_DATA);
    % Replace missing data value with NaN
    data{1}(data{1} == MIS_DATA) = nan;

    B = data{1}(:,1:3);
    E = data{1}(:,4:5);

    B(B == MIS_DATA) = nan;
    
    E(abs(E)>30) = nan; % Obvious spikes

    logmsg('Saving %s\n',matfile);
    save(matfile,'E','B','Header');
else
    logmsg(sprintf('Reading %s\n',matfile));
    load(matfile);
end
