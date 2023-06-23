%% Z nonzero at single frequency on DFT grid (so no leakage)
%  with time domain window

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
k  = 25;
Z  = (1+1j)/sqrt(2);
f  = k/Nt;
wf = 1;

regstr = sprintf('OLS/$N_b=%d$',2*wf+1); 

Sx_opts = struct('Nt',Nt,'Z',Z,'f',f,'dB',0.0,'dE',0.0);
Sx = demo_signals('simple',Sx_opts);
Sx.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.window.function = @tdwindow;
    opts1.td.window.functionstr = 'Rectangle';
    opts1.td.window.functionargs = {@rectwin};
S1 = tflab(Sx.In,Sx.Out,opts1);
S1.Options.description = sprintf('rectwin/%s',regstr);

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.window.function = @tdwindow;
    opts2.td.window.functionstr = 'Parzen';
    opts2.td.window.functionargs = {@parzenwin};
S2 = tflab(Sx.In,Sx.Out,opts2);
S2.Options.description = sprintf('parzenwin/%s',regstr);

dockreset();

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','windowed'));
figure();
    tsplot(S2,struct('type','windowed'));
figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure();
    dftplot({S1,S2},struct('type','final-raw'));
figure();
    dftplot({S1,S2},struct('type','error-raw'));

figure();
    zplot({Sx,S1,S2});
figure();
    zplot({Sx,S1,S2},struct('type',3));

if wf > 0    
    figure();
        qqplot_(S1,k);
end
