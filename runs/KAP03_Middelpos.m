close all;
set(0,'defaultFigureWindowStyle','docked');

lemi  = 0; 
savefmt = struct('pdf',0,'png',1);

if 1
    cid   = 'Middelpos';
    id    = 'Middelpos';
    start = '2012-07-12T00:00:00.000'; % Must be in this format
    stop  = '2012-09-04T00:00:00.000'; % Must be in this format
    timedelta = 1;
    % Read input/output data
    [B,E] = Middelpos_data();
end

if 1
    cid   = 'KAP03';
    id    = 'KAP103';
    %id    = 'KAP106';
    %id    = 'KAP109';
    
    % Read input/output data
    [B,E,H] = KAP03_data(id);
    start = [strrep(H.STARTTIME,' ','T'),'.000'];
    stop  =  [strrep(H.ENDTIME,' ','T'),'.000'];
    timedelta = str2num(H.DELTA_T);
end

addpath([fileparts(mfilename('fullpath')),'/../window']); 
addpath([fileparts(mfilename('fullpath')),'/../fft']); 
addpath([fileparts(mfilename('fullpath')),'/../lib']); 
addpath([fileparts(mfilename('fullpath')),'/../']);
addpath([fileparts(mfilename('fullpath')),'/../plot']);
    
% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

% Prefix for figure files
figprfx = [scriptpath,sprintf('/figures/%s/%s-',cid,id)];

% Ouput file
matfile = [scriptpath,sprintf('/data/%s/%s',cid,id)];

% Plot title information
ptitle = sprintf('%s %s - %s',id,start,stop);

% Variable name information
iopts = struct('info',struct(),'td',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/m';
iopts.info.timeunit = 's';
iopts.info.timedelta = timedelta;
iopts.info.timestart = start;

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
    S3 = transferfnFD_lemimt(B,E);
    S3.In = B(:,1:2);
    S3.Out = E;
end

% Other code only uses Bx and By.
B = B(:,1:2);

if 1
    desc1 = sprintf('%d 1-day segments',size(B,1)/ppd);
    opts1 = transferfnFD_options(1,iopts);
        opts1.transferfnFD.loglevel = 1;
        opts1.td.window.width = 86400/timedelta;
        opts1.td.window.shift = 86400/timedelta;
end

if 1
    desc2 = sprintf('One %d-day segment',size(B,1)/ppd);
    opts2 = transferfnFD_options(1,iopts);
        opts2.transferfnFD.loglevel = 1;
end

% Execute runs
%I = 1:86400*10; % Use a sub-set of data
I = 1:size(B,1); % Use all 55 days of data
S1 = transferfnFD(B(I,:),E(I,:),opts1);
S2 = transferfnFD(B(I,:),E(I,:),opts2);

% Test S2.Z on same segments as S1.
S2 = transferfnMetrics(S2,opts2,S1.Segment.IndexRange);

if lemi
    % Test S3.Z on same segments as S1.
    S3 = transferfnMetrics(S3,opts2,S1.Segment.IndexRange);
    S3.Options = struct('description','LEMI One 55-day segment');
end

% Modify default descriptions of run
S1.Options.description = desc1;
S2.Options.description = desc2;

% Restrict period range in plots
band = [f(Ic(2)),f(Ic(end))];
period_range = [1/band(2),1/band(1)]*iopts.info.timedelta;


timeseries_plot(S1,...
    struct('title',ptitle,'filename',[figprfx,'timeseries_1'],...
           'savefmt',savefmt));
timeseries_plot(S1,...
    struct('type','error','title',[ptitle,'; ',desc1],...
           'filename',[figprfx,'timeseries-error_1'],'savefmt',savefmt));
timeseries_plot(S2,...
    struct('type','error','title',[ptitle,'; ',desc1],...
           'filename',[figprfx,'timeseries-error_2'],'savefmt',savefmt));

spectrum_plot(S1,...
    struct('type','raw','title',[ptitle,'; ',desc1],...
           'period_range',period_range,...
           'filename',[figprfx,'spectrum_1'],'savefmt',savefmt));
spectrum_plot(S2,...
    struct('type','raw','title',[ptitle,'; ',desc1],...
           'period_range',period_range,...
           'filename',[figprfx,'spectrum_2'],'savefmt',savefmt));

sn_plot(S1,...
    struct('title',[ptitle,'; ',desc1],'period_range',period_range,...
           'filename',[figprfx,'SN_1'],'savefmt',savefmt));
sn_plot(S2,...
    struct('title',[ptitle,'; ',desc1],'period_range',period_range,...
           'filename',[figprfx,'SN_2'],'savefmt',savefmt));

filename = [figprfx,'SN'];
if lemi    
    sn_plot({S1,S2,S3},...
        struct('period_range',period_range,'title',ptitle,...
               'filename',filename,'savefmt',savefmt));
else
    sn_plot({S1,S2},...
        struct('period_range',period_range,'title',ptitle,...
               'filename',filename,'savefmt',savefmt));
end

transferfnZ_plot(S1,...
    struct('period_range',period_range,'title',ptitle,...
    'filename',[figprfx,'Z_1'],'savefmt',savefmt));
transferfnZ_plot(S1,...
    struct('period_range',period_range,'title',ptitle,...
    'filename',[figprfx,'Z_2'],'savefmt',savefmt));

filename = [figprfx,'Z'];    
if lemi    
    transferfnZ_plot({S1,S2,S3},...
        struct('period_range',period_range,'title',ptitle,...
               'filename',filename,'savefmt',savefmt));
else
    transferfnZ_plot({S1,S2},...
        struct('period_range',period_range,'title',ptitle,...
               'filename',filename,'savefmt',savefmt));
end

if lemi
    S = {S1,S2,S3};
else
    S = {S1,S2};
end
save(matfile,'-v7.3','S');
