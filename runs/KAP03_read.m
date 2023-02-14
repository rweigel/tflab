function [B,E,t,matfile,Header] = KAP03_read(id)

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

% Dir of data file
datapath = fullfile(scriptdir(),'data','KAP03'); 

% Content of data file
matfile = fullfile(datapath,'/',[id,'_raw.mat']);

if exist(matfile,'file')
    logmsg(sprintf('Reading %s\n',matfile));
    load(matfile);
    %return
end

fname = fullfile(datapath,'measurements',[lower(id),'as.ts']);
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

start = [strrep(Header.STARTTIME,' ','T'),'.000'];
dno = datenum(datevec(start,'yyyy-mm-ddTHH:MM:SS.FFF'));
t = dno + 5*[1:size(E,1)]/86400;

% Verify that time series is on uniform grid
assert(length(t) == size(B,1));
assert(length(t) == size(E,1));

logmsg('Saving %s\n',matfile);
save(matfile,'E','B','t','Header');
