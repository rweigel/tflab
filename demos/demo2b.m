%% Single frequency not on DFT grid (so leakage)
% Demonstrates use of non-rectangular time domain window (Parzen).

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
f  = 25.5/Nt;
wf = 0;

regstr = sprintf('OLS/$N_b=%d$',2*wf+1); 

Sx_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sx = demo_signals('simple',Sx_opts);
Sx.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.window.function = @tdwindow;
    opts1.td.window.functionstr = 'rectwin';
    opts1.td.window.functionargs = {@rectwin};
S1 = tflab(Sx.In,Sx.Out,opts1);
S1.Options.description = sprintf('rectwin/%s',regstr);

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.window.function = @tdwindow;
    opts2.td.window.functionstr = 'parzenwin';
    opts2.td.window.functionargs = {@parzenwin};
    %opts2.td.window.functionstr = 'DPSS; $N_w=2$, $K=1$';
    %opts2.td.window.functionargs = {@dpss, 2, 1};
    
S2 = tflab(Sx.In,Sx.Out,opts2);
S2.Options.description = sprintf('%s/%s',regstr,opts2.td.window.functionstr);

dockreset();

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S2,struct('type','windowed'));
figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure();
    dftplot({S1,S2},struct('type','windowed'));
figure();
    dftplot({S1,S2},struct('type','error'));

figure();
    snplot(S1);
    
figure();
    zplot({Sx,S1,S2});
figure();
    zplot({Sx,S1,S2},struct('type',3));

if wf > 0
    figure();
        qqplot_(S1,Nf-wf+1);
end
