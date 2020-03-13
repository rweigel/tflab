if ~exist('Middelpos.mat')

    dname = '/home/weigel/git/lemimt/data/2012/all';
    do = datenum(2012,7,12);
    df = datenum(2012,9,4);

    D = [];
    for d = do:df
        [y,m,d] = datevec(d);
        fname =sprintf('%s/MT_Middelpos_%d%02d%02d.t82',dname,y,m,d);
        fprintf('Reading %s\n',fname);
        tmp = load(fname);
        D = [D;tmp];
    end
    B = D(:,7:9);
    E = D(:,12:13);
    t = datenum(D(:,1:6));
    save Middelpos.mat D B E t
else
    load Middelpos
end

pdf = 1;

B = B(:,1:2);
[B,E] = removemean(B,E);

% Plot title information
ptitle = 'Middelpos 2012-07-12 - 2012-09-04';

iopts = struct('info',struct());
iopts.info.instr = {'$B_x$','$B_y$'};
iopts.info.inunit= 'nT';
iopts.info.outstr = {'$E_x$','$E_y$'};
iopts.info.outunit= 'mV/m';
iopts.info.timeunit = 's';
iopts.info.timestart = '2012-07-12T00:00:00.000';

if 1
desc1 = 'Ave. of 55 1-day segments';
opts1 = transferfnFD_options(1,iopts);
    opts1.transferfnFD.loglevel = 1;
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
    
S1 = transferfnFD(B,E,opts1);
S2 = transferfnFD(B,E,opts2);

S1.Options.description = desc1;
S2.Options.description = desc2;

close all;
figure();clf;
    timeseries_plot(S1,struct('title',ptitle));
    %figsave(pdf,'Middelpos-timeseries_1.pdf');
    
figure();clf;
    timeseries_plot(S1,struct('title',[ptitle,'Parzen tapering'],'type','windowed'));
    %figsave(pdf,'Middelpos-timeseries_windowed_1.pdf');
    
figure();clf;
    spectrum_plot(S1,struct('type','raw','title',[ptitle,'; ',desc1]));
    %figsave(pdf,'Middelpos-spectrum_1.pdf');        

figure();clf;
    spectrum_plot(S1,struct('type','windowed','title',[ptitle,'; ',desc1]));
    %figsave(pdf,'Middelpos-spectrum_windowed_1.pdf');

figure();clf;
    spectrum_plot(S2,struct('type','raw','title',[ptitle,'; ',desc2]));
    %figsave(pdf,'Middelpos-spectrum_2.pdf');

figure(4);clf;
    timeseries_plot(S1,struct('type','error','title',[ptitle,'; ',desc1]));
    %figsave(pdf,'Middelpos-error_1.pdf');

figure(5);clf;
    timeseries_plot(S2,struct('type','error','title',[ptitle,'; ',desc2]));
    %figsave(pdf,'Middelpos-error_2.pdf');  
    
period_range = [min(1./S1.fe),max(1./S1.fe(2:end))];

figure(6);clf;
    sn_plot(S1,S2,struct('period_range',period_range,'title',ptitle));
    figsave(pdf,'Middelpos-error.pdf');

figure(7);clf;
    popts.filename = 'Middelpos-Z';
    transferfnZ_plot(S1,S2,struct('period_range',period_range,'title',ptitle));
