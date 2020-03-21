pdf = 0; % Save plots to pdf file.
png = 0; % Save plots to png file.

% Dir of this script
scriptpath = fileparts(mfilename('fullpath')); 

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

% Read input/output data
[B,E,t] = Middelpos_data();
[B,E] = removemean(B,E);

% Keep frequencies in range of minimum to maximum frequency used to
% estimate TF when 1-day segments are used.
[~,f] = fftfreq(86400);
[fe,Ic,Ne] = evalfreq(86400);
band = [f(Ic(2)-Ne(2)),f(Ic(end)+Ne(end))];
B = bandpass(B,band);
E = bandpass(E,band);

if 1
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
    desc1 = '55 1-day segments';
    opts1 = transferfnFD_options(1,iopts);
        opts1.transferfnFD.loglevel = 1;
        opts1.td.window.width = 86400;
        opts1.td.window.shift = 86400;
end

if 1
    desc2 = 'One 55-day segment';
    opts2 = transferfnFD_options(1,iopts);
        opts2.transferfnFD.loglevel = 1;
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

% Execute runs
%I = 1:86400*10; % Use a sub-set of data from 55 days.
I = 1:size(B,1); % Use all 55 days of data
S1 = transferfnFD(B(I,:),E(I,:),opts1);
S2 = transferfnFD(B(I,:),E(I,:),opts2);

% Test S2.Z on same segments as S1.
S2 = transferfnMetrics(S2,opts2,S1.Segment.IndexRange);

if 1
    % Test S3.Z on same segments as S1.
    S3 = transferfnMetrics(S3,opts2,S1.Segment.IndexRange);
    S3.Options = struct('description','LEMI One 55-day segment');
end

% Modify default descriptions of run
S1.Options.description = desc1;
S2.Options.description = desc2;

% Restrict period range in plots
period_range = [1/band(2),1/band(1)];

close all;
filename = [scriptpath,'/figures/Middelpos/Middelpos-SN.pdf'];
timeseries_plot(S1,struct('title',ptitle,'filename',filename));

filename = [scriptpath,'/figures/Middelpos/Middelpos-SN.pdf'];
    timeseries_plot(S1,struct('type','error','title',[ptitle,'; ',desc1]));
    figsave(pdf,'/figures/Middelpos/Middelpos-timeseries-error_1.pdf');

figure();clf;
    timeseries_plot(S2,struct('type','error','title',[ptitle,'; ',desc2]));
    figsave(pdf,'/figures/Middelpos/Middelpos-timeseries-error_2.pdf');  
        
figure();clf;
    spectrum_plot(S1,struct('type','raw','title',[ptitle,'; ',desc1],'period_range',period_range));
    figsave(pdf,'/figures/Middelpos/Middelpos-spectrum_1.pdf');        

figure();clf;
    spectrum_plot(S2,struct('type','raw','title',[ptitle,'; ',desc2],'period_range',period_range));
    figsave(pdf,'/figures/Middelpos/Middelpos-spectrum_2.pdf');

figure();clf;
    filename = [scriptpath,'/figures/Middelpos/Middelpos-SN_1.pdf'];
    sn_plot(S1,struct('period_range',period_range,'title',[ptitle,'; ',desc1],'filename',filename));

figure();clf;
    filename = [scriptpath,'/figures/Middelpos/Middelpos-SN_2.pdf'];
    sn_plot(S2,struct('period_range',period_range,'title',[ptitle,'; ',desc2],'filename',filename));
    
filename = [scriptpath,'/figures/Middelpos/Middelpos-SN.pdf'];
sn_plot({S1,S2,S3},struct('period_range',period_range,'title',ptitle,'filename',filename));
    
filename = [scriptpath,'/figures/Middelpos/Middelpos-Z_1.pdf'];
transferfnZ_plot(S1,struct('period_range',period_range,'title',ptitle,'filename',filename));

filename = [scriptpath,'/figures/Middelpos/Middelpos-Z_2.pdf'];
transferfnZ_plot(S1,struct('period_range',period_range,'title',ptitle,'filename',filename));
    
filename = [scriptpath,'/figures/Middelpos/Middelpos-Z.pdf'];
transferfnZ_plot({S1,S2,S3},struct('period_range',period_range,'title',ptitle,'filename',filename));
