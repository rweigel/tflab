function [B,E,t,infile,outfile,timedelta] = KAP03_clean(sta)

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();
addpath(fullfile(fileparts(mfilename('fullpath')),'readers'));

% Dir of data file
outfile = fullfile(scriptdir(),'data','KAP03',sta,[sta,'_cleaned.mat']); 

if exist(outfile,'file')
    logmsg('Reading %s\n',outfile);
    load(outfile)
    %return
end

% Read input/output data
[B,E,t,infile,Header] = KAP03_read(sta);

timedelta = str2num(Header.DELTA_T);

start = [strrep(Header.STARTTIME,' ','T'),'.000'];
%stop  = [strrep(Header.ENDTIME,' ','T'),'.000'];
dno = datenum(datevec(start,'yyyy-mm-ddTHH:MM:SS.FFF'));

figprep();
subplot(3,1,1)
    plot(t,E)
    title('Raw');
    datetick();
    grid on
    set(gca,'XTickLabels',[])
subplot(3,1,2)
    E = despikeE(E, sta);
    plot(t,E)
    title('Despiked');
    datetick()
    grid on
    set(gca,'XTickLabels',[])    
subplot(3,1,3)
    E = naninterp1(E);
    plot(t,E)
    title('Interpolated');
    datetick()
    grid on
    xlabel('Month/Day of 2003')
figsave([outfile(1:end-3),'png']);

logmsg('Interpolating over NaNs in B\n');
B = naninterp1(B);

function X = despikeE(X, sta)
    Ibad = [];
    if strcmp(sta,'KAP103')
        % From visual inspection
        Ibad = [242900, 242910;...
                308981, 309282;...
                341328, 341950;...
                419932, 422771];
    end
    if strcmp(sta,'KAP163')
        Ibad = [188720, 188760];
    end
    for i = 1:size(Ibad,1)
        X(Ibad(i,1):Ibad(i,2),:) = nan;
    end
end

logmsg('Saving %s\n',outfile);
save(outfile,'B','E','t','infile','outfile','timedelta','Header')

end

