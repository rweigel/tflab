function Middelpos_main(rundir, filestr, start, stop, Nboot)

% Read input/output data
[B,E,t] = Middelpos_clean(start,stop,rundir,0);

% Make length an integer number of segments.
pps = 86400;
I = pps*floor(size(B,1)/pps);
B = B(1:I,:);
E = E(1:I,:);
t = t(1:I);

%% Set common metadata
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = 1;
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.timestop  = datestr(t(end),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.chainid   = 'SANSA';
    meta.stationid = 'Middelpos';

tfn = 0;
for const_term = 0:1
    tfn = tfn+1;
    logmsg('-- Computing TF%d --\n',tfn);
    ppd = 86400;
    if const_term == 0
        desc = sprintf('1 %d-day segment',size(B,1)/ppd);
    else
        desc = sprintf('1 %d-day segment; $\\delta$ term',size(B,1)/ppd);
    end
    opts{tfn} = tflab_options(1);
        opts{tfn}.tflab.loglevel = 1;
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    if const_term == 0
        [~,t_latex] = evalfreq_log(size(B,1), opts{tfn}.fd.evalfreq.functionargs{:});
        fid = fopen(fullfile(rundir,append(opts{tfn}.filestr,'-evalfreq_table.tex')),'w');
        fprintf(fid, t_latex);
        fclose(fid);
    end
    
    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
end

for const_term = 0:1
    tfn = tfn+1;
    logmsg('-- Computing TF%d --\n',tfn);
    pps = 86400;
    if const_term == 0
        desc = sprintf('%d %d-day segments',size(B,1)/pps,pps/86400);
    else
        desc = sprintf('%d %d-day segments; $\\delta$ term',size(B,1)/pps,pps/86400);
    end
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.window.width = pps;
        opts{tfn}.td.window.shift = pps;
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    if const_term == 0        
        [~,t_latex] = evalfreq_log(size(B,1), opts{tfn}.fd.evalfreq.functionargs{:});
        fid = fopen(fullfile(rundir,append(opts{tfn}.filestr,'-evalfreq_table.tex')),'w');
        fprintf(fid, t_latex);
        fclose(fid);
    end

    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
end

if 1
    const_term = 0;

    tfn = tfn+1;
    logmsg('-- Computing TF%d --\n',tfn);
    Tm = 1*86400;
    band = [1/Tm,0.5];
    desc = '1 %d-day seg.; %d-day bandpass;';
    desc = sprintf(desc,size(B,1)/pps,Tm/86400);
    
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.detrend.function = @bandpass_;
        opts{tfn}.td.detrend.functionstr = 'bandpass_';
        opts{tfn}.td.detrend.functionargs = {band};
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
    
    tfn = tfn+1;
    logmsg('-- Computing TF%d --\n',tfn);
    pps = 86400;
    Tm = 1*86400;
    band = [1/Tm,0.5];
    desc = '%d %d-day segments; %d-day bandpass';
    desc = sprintf(desc,size(B,1)/pps,pps/86400,Tm/86400);
    
    opts{tfn} = tflab_options(1);
        opts{tfn}.td.window.width = pps;
        opts{tfn}.td.window.shift = pps;
        opts{tfn}.td.detrend.function = @bandpass_;
        opts{tfn}.td.detrend.functionstr = 'bandpass_';
        opts{tfn}.td.detrend.functionargs = {band};
        opts{tfn}.fd.regression.const_term = const_term;
        opts{tfn}.description = desc;
        opts{tfn}.fd.bootstrap.N = Nboot;
        opts{tfn}.fd.bootstrap.nmin = size(B,1)/pps;
        opts{tfn}.filestr = sprintf('%s-tf%d',filestr,tfn);
    
    TFs{tfn} = tflab(B(:,1:2),E,opts{tfn});
    TFs{tfn}.Metadata = meta;
    
    savetf(TFs{tfn}, fullfile(rundir,opts{tfn}.filestr));
    TFs{tfn} = [];
end
