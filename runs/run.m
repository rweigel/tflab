function S = run(sid,lemi)

zread_dir = [fileparts(mfilename('fullpath')),'/zread'];
if ~exist(zread_dir,'dir')
    url = 'https://github.com/rweigel/zread';
    com = sprintf('cd %s; git clone %s; cd zread; git checkout 8be764e50439db308bbb0b51b886bf0b7fb10c24',fileparts(mfilename('fullpath')),url);
    fprintf('Calling system with command %s\n',com);
    [status,msg] = system(com);
    if status ~= 0
        error('System command failed: %s\nMessage:\n%s\n',com,msg);
    end
end
addpath(zread_dir);

if strcmp(sid,'Middelpos')
    cid = 'Middelpos';  % Chain ID
    sid = 'Middelpos'; % Station ID

    [B,E] = data_Middelpos(); % Read input/output data

    start = '2012-07-12T00:00:00.000'; % Must be in this format
    stop  = '2012-09-04T00:00:00.000'; % Must be in this format
    timedelta = 1;
end

if strncmp(sid,'KAP',3)
    cid = 'KAP03';
    
    % Read input/output data
    [B,E,Info] = data_KAP03(sid);
    start = [strrep(Info.STARTTIME,' ','T'),'.000'];
    stop  =  [strrep(Info.ENDTIME,' ','T'),'.000'];
    timedelta = str2num(Info.DELTA_T);
end
%% End configuration

%% Run
addpath([fileparts(mfilename('fullpath')),'/../']);
addpath([fileparts(mfilename('fullpath')),'/../window/']);
addpath([fileparts(mfilename('fullpath')),'/../fft/']);
addpath([fileparts(mfilename('fullpath')),'/../lib/']);

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Ouput file
matfile = [scriptpath,sprintf('/data/%s/%s',cid,sid)];

% Plot title information
ptitle = sprintf('%s %s - %s',sid,start,stop);

% Variable name information
iopts = struct('info',struct(),'td',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/m';
iopts.info.timeunit = 's';
iopts.info.timedelta = timedelta;
iopts.info.timestart = start;
iopts.info.stationid = sid;
iopts.info.chainid = cid;

%% Remove mean, interpolate over NaNs, pad E if NaNs at start or end
[B,E] = removemean(B,E);

% TODO: Report on largest gap
ti = [1:size(B,1)]';
for i = 1:size(B,2)
    tg = find(~isnan(B(:,i)));
    B(:,i) = interp1(tg,B(tg,i),ti);
end
for i = 1:size(E,2)
    tg = find(~isnan(E(:,i)));    
    E(:,i) = interp1(tg,E(tg,i),ti);    
end

for i = 1:size(E,2)
    I = find(isnan(E(:,i)));
    if length(I) > 0
        logmsg(sprintf('Set %d leading or trailing NaNs in column %d of E to zero\n',length(I),i))
        E(I,:) = 0;
    end
end

%% Trim so integer number of days
% TODO: This assumes segment length will be one day. Generalize.
ppd = 86400/iopts.info.timedelta;

Ix = floor(size(B,1)/ppd);
B = B(1:ppd*Ix,:);
E = E(1:ppd*Ix,:);

if 0
    %% Band pass
    % Keep frequencies in range of minimum to maximum frequency used to
    % estimate TF when 1-day segments are used.
    [~,f] = fftfreq(ppd);
    [fe,Ic,Ne] = evalfreq(ppd);
    band = [f(Ic(2)-Ne(2)),f(Ic(end)+Ne(end))];
    B = bandpass(B,band);
    E = bandpass(E,band);
end

%% Compute TFs

% First TF
desc1 = sprintf('OLS; %d 1-day segments',size(B,1)/ppd);
opts1 = transferfnFD_options(1,iopts);
    opts1.transferfnFD.loglevel = 1;
    opts1.td.window.width = 86400/timedelta;
    opts1.td.window.shift = 86400/timedelta;

S{1} = transferfnFD(B(:,1:2),E,opts1);
% Modify default descriptions of run
S{1}.Options.description = desc1;

%% Second TF
desc2 = sprintf('OLS; One %d-day segment',size(B,1)/ppd);
opts2 = transferfnFD_options(1,iopts);
    opts2.transferfnFD.loglevel = 1;

S{2} = transferfnFD(B(:,1:2),E,opts2);           
S{2}.Options.description = desc2;
    
% Test S{2}.Z on same segments as S{1}.
S{2} = transferfnMetrics(S{2},opts2,S{1}.Segment.IndexRange);

tf = length(S);

%% TF computed using BIRP
if strcmp(cid,'KAP03')
    tf = tf+1;
    S{tf} = read_edi([zread_dir,'/data/kap003.edi']);
    S{tf}.Z(S{tf}.Z > 1e31) = NaN;
    
    S{tf}.Options.info = S{1}.Options.info;
    S{tf}.Options.description = 'BIRP';
    S{tf}.In  = S{1}.In;
    S{tf}.Out = S{1}.Out;
    S{tf}.Time = S{1}.Time;
    S{tf} = transferfnMetrics(S{tf},S{1}.Options,S{1}.Segment.IndexRange,1);
end

%% Compute TF using LEMI MT
if lemi
    tf = tf+1;
    % Use LEMI MT program to compute TF. Note that LEMI MT requires
    % three components of B. (Bz is used to compute another TF that is
    % not used here).
    addpath(lemimt_dir);
    S{tf} = transferfnFD_lemimt(B,E);
    S{tf}.In = B(:,1:2);
    S{tf}.Out = E;
    % Test S{tf}.Z on same segments as S{1}.
    S{tf}.Time = S{1}.Time;
    S{tf} = transferfnMetrics(S{tf},opts2,S{1}.Segment.IndexRange);
    S{tf}.Options.info = S{1}.Options.info;
    S{tf}.Options.description = 'LEMI; One 55-day segment';
end

save(matfile,'-v7.3','S');
