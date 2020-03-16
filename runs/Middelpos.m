pdf = 0; % Save plots to pdf file.
png = 0; % Save plots to png file.

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

[B,E] = removemean(B,E);

fe = evalfreq(86400);
B = bandpass(B,[fe(2),fe(end)]);
E = bandpass(E,[fe(2),fe(end)]);

if 1
    addpath('/home/weigel/git/lemimt');
    S3 = transferfnFD_lemimt(B,E);
    S3.In = B(:,1:2);
    S3.Out = E;
end

B = B(:,1:2);
% Plot title information
ptitle = 'Middelpos 2012-07-12 - 2012-09-04';

% Variable name information
iopts = struct('info',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/m';
iopts.info.timeunit = 's';
iopts.info.timestart = '2012-07-12T00:00:00.000';

if 1
    desc1 = '55 1-day segments';
    opts1 = transferfnFD_options(1,iopts);
        opts1.transferfnFD.loglevel = 1;
        opts1.td.window.width = 86400;
        opts1.td.window.shift = 86400;
end

if 0
    desc1 = '55 1-day segments; No segment averaging.';
    opts1 = transferfnFD_options(1,iopts);
        opts1.transferfnFD.loglevel = 1;
        opts1.fd.stack.average.function = '';
        opts1.td.window.width = 86400;
        opts1.td.window.shift = 86400;
end

if 0
    desc1 = 'One 55-day segment; Parzen tapering.';
    opts1 = transferfnFD_options(1,iopts);
        opts1.td.window.function = @tdwindow; 
        opts1.td.window.functionstr = 'Parzen';        
        opts1.td.window.functionargs = {@parzenwin};        
        opts1.transferfnFD.loglevel = 1;
end

desc2 = 'One 55-day segment';
opts2 = transferfnFD_options(1,iopts);
    opts2.transferfnFD.loglevel = 1;

% Execute runs
%I = 1:86400*10;
I = 1:size(B,1);
S1 = transferfnFD(B(I,:),E(I,:),opts1);
S2 = transferfnFD(B(I,:),E(I,:),opts2);

% Test S2.Z on same segments as S1.
S2 = transferfnMetrics(S2,opts2,S1.Segment.IndexRange);

if 1
    S3.In = B(:,1:2);
    S3.Out = E;
    S3 = transferfnMetrics(S3,opts2,S1.Segment.IndexRange);
    S3.Options = struct('description','LEMI One 55-day segment');
end

% Modify default descriptions of run
S1.Options.description = desc1;
S2.Options.description = desc2;

period_range = [min(1./S1.fe),max(1./S1.fe(2:end))];

close all;
figure();clf;
    timeseries_plot(S1,struct('title',ptitle));
    figsave(pdf,'Middelpos-timeseries_1.pdf');
    
figure();clf;
    timeseries_plot(S1,struct('type','error','title',[ptitle,'; ',desc1]));
    figsave(pdf,'Middelpos-timeseries-error_1.pdf');

figure();clf;
    timeseries_plot(S2,struct('type','error','title',[ptitle,'; ',desc2]));
    figsave(pdf,'Middelpos-timeseries-error_2.pdf');  
    
%figure();clf;
    %timeseries_plot(S1,struct('title',[ptitle,'Parzen tapering'],'type','windowed'));
    %figsave(pdf,'Middelpos-timeseries_windowed_1.pdf');
    
figure();clf;
    spectrum_plot(S1,struct('type','raw','title',[ptitle,'; ',desc1],'period_range',period_range));
    figsave(pdf,'Middelpos-spectrum_1.pdf');        

%figure();clf;
    %spectrum_plot(S1,struct('type','windowed','title',[ptitle,'; ',desc1]));
    %figsave(pdf,'Middelpos-spectrum_windowed_1.pdf');

figure();clf;
    spectrum_plot(S2,struct('type','raw','title',[ptitle,'; ',desc2],'period_range',period_range));
    figsave(pdf,'Middelpos-spectrum_2.pdf');

figure();clf;
    filename = [scriptpath,'/figures/Middelpos/Middelpos-SN.pdf'];
    sn_plot({S1,S2,S3},struct('period_range',period_range,'title',ptitle,'filename',filename));

figure();clf;
    filename = [scriptpath,'/figures/Middelpos/Middelpos-Z_1.pdf'];
    transferfnZ_plot(S1,struct('period_range',period_range,'title',ptitle,'filename',filename));
    
figure();clf;
    filename = [scriptpath,'/figures/Middelpos/Middelpos-Z.pdf'];
    transferfnZ_plot({S1,S2,S3},struct('period_range',period_range,'title',ptitle,'filename',filename));
