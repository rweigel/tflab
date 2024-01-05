function [B,E,t,infile,outfile,timedelta] = KAP03_clean(sta,plot)

if nargin < 2
    plot = 0;
end

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

if plot
    Er = E;
    Ed = despikeE(Er, sta);
else
    Ed = despikeE(E, sta);
end
logmsg('Interpolating over NaNs in E\n');
E = naninterp1(Ed);

logmsg('Interpolating over NaNs in B\n');
B = naninterp1(B);

if plot
    figprep();
    subplot(3,1,1)
        plot(t,Er)
        title('Raw');
        datetick();
        grid on
        set(gca,'XTickLabels',[])
    subplot(3,1,2)
        plot(t,Ed)
        title('Despiked');
        datetick()
        grid on
        set(gca,'XTickLabels',[])    
    subplot(3,1,3)
        plot(t,E)
        title('Interpolated');
        datetick()
        grid on
        xlabel('Month/Day of 2003')
    figsave([outfile(1:end-3),'png']);
end

logmsg('Saving %s\n',outfile);
save(outfile,'B','E','t','infile','outfile','timedelta','Header')

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
        Ibad = [188720, 188760; 152629, 152708];
    end
    for i = 1:size(Ibad,1)
        X(Ibad(i,1):Ibad(i,2),:) = nan;
    end
end


end

