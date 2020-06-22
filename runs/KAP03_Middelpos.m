clear

lemi  = 0;

if 0
    cid   = 'Middelpos';
    sid    = 'Middelpos';
    start = '2012-07-12T00:00:00.000'; % Must be in this format
    stop  = '2012-09-04T00:00:00.000'; % Must be in this format
    timedelta = 1;
    % Read input/output data
    [B,E] = Middelpos_data();
end

if 1
    cid   = 'KAP03'; 
    sid    = 'KAP103';
    %id    = 'KAP106';
    %id    = 'KAP109';
    
    % Read input/output data
    [B,E,H] = KAP03_data(sid);
    start = [strrep(H.STARTTIME,' ','T'),'.000'];
    stop  =  [strrep(H.ENDTIME,' ','T'),'.000'];
    timedelta = str2num(H.DELTA_T);
end

addpath([fileparts(mfilename('fullpath')),'/../']);
addpath([fileparts(mfilename('fullpath')),'/../window']); 
addpath([fileparts(mfilename('fullpath')),'/../fft']); 
addpath([fileparts(mfilename('fullpath')),'/../lib']); 
addpath([fileparts(mfilename('fullpath')),'/../plot']);
    
% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Prefix for figure files
figprfx = [scriptpath,sprintf('/figures/%s/%s-',cid,sid)];

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

ppd = 86400/iopts.info.timedelta;

Ix = floor(size(B,1)/ppd);
% Trim so integer number of segments
B = B(1:ppd*Ix,:);
E = E(1:ppd*Ix,:);

% Keep frequencies in range of minimum to maximum frequency used to
% estimate TF when 1-day segments are used.
[~,f] = fftfreq(ppd);
[fe,Ic,Ne] = evalfreq(ppd);
band = [f(Ic(2)-Ne(2)),f(Ic(end)+Ne(end))];
B = bandpass(B,band);
E = bandpass(E,band);

if lemi
    % Use LEMI MT program to compute TF. Note that LEMI MT requires
    % three components of B. (Bz is used to compute another TF that is
    % not used here).
    addpath('/home/weigel/git/lemimt');
    S{3} = transferfnFD_lemimt(B,E);
    S{3}.In = B(:,1:2);
    S{3}.Out = E;
end

% Other code only uses Bx and By.
B = B(:,1:2);

desc1 = sprintf('OLS; %d 1-day segments',size(B,1)/ppd);
opts1 = transferfnFD_options(1,iopts);
    opts1.transferfnFD.loglevel = 1;
    opts1.td.window.width = 86400/timedelta;
    opts1.td.window.shift = 86400/timedelta;

desc2 = sprintf('OLS; One %d-day segment',size(B,1)/ppd);
opts2 = transferfnFD_options(1,iopts);
    opts2.transferfnFD.loglevel = 1;

% Execute runs
%I = 1:86400*10; % Use a sub-set of data
I = 1:size(B,1); % Use all 55 days of data
S{1} = transferfnFD(B(I,:),E(I,:),opts1);
S{2} = transferfnFD(B(I,:),E(I,:),opts2);           

% Modify default descriptions of run
S{1}.Options.description = desc1;
S{2}.Options.description = desc2;

% Test S{2}.Z on same segments as S{1}.
S{2} = transferfnMetrics(S{2},opts2,S{1}.Segment.IndexRange);

if lemi
    % Test S{3}.Z on same segments as S{1}.
    S{3}.Time = S{1}.Time;
    S{3} = transferfnMetrics(S{3},opts2,S{1}.Segment.IndexRange);
    S{3}.Options = struct('description','LEMI; One 55-day segment');
end

if strcmp(cid,'KAP03')
    zread = '/Users/robertweigel/git/zread/';
    addpath(zread);
    S{4} = read_edi([zread,'data/samtex/Final_with_DeBeers/edi_files/KAP03/kap003.edi']);
    S{4}.Z(S{4}.Z > 1e31) = NaN;
    
    % Does not matter
    %I = find(S{4}.fe < 0.5);
    %S{4}.fe = S{4}.fe(I);    
    %S{4}.Z = S{4}.Z(I,:);

    S{4}.Options.info = S{1}.Options.info;
    S{4}.Options.description = 'BIRP';
    S{4}.In  = S{1}.In;
    S{4}.Out = S{1}.Out;
    S{4}.Time = S{1}.Time;
    S{4} = transferfnMetrics(S{4},S{1}.Options,S{1}.Segment.IndexRange,1);
end

save(matfile,'-v7.3','S');

