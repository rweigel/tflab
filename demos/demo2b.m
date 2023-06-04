%% Single frequency not on DFT grid (so leakage)
% Demonstrates use of non-rectangular time domain window (Parzen).

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
f  = 25.5/Nt;
wf = 0;

Sx_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sx = demo_signals('simple',Sx_opts);

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.window.function = @tdwindow;
    opts1.td.window.functionstr = 'Rectangle';
    opts1.td.window.functionargs = {@rectwin};
S1 = tflab(Sx.In,Sx.Out,opts1);
S1.Options.description = 'Rectangle';

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.window.function = @tdwindow;
    opts2.td.window.functionstr = 'Parzen';
    opts2.td.window.functionargs = {@dpss, 1, 1};
S2 = tflab(Sx.In,Sx.Out,opts2);
S2.Options.description = 'Parzen';

% Use same variable labels from S1 for Sa
Sx.Options.info = S1.Options.info; 
Sx.Options.description = 'Actual';


dock on;figure(1);close all;

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','windowed'));
figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure();
    dftplot({S1,S2},struct('type','windowed'));
figure();
    dftplot({S1,S2},struct('type','error'));

figure();
    zplot({Sx,S1,S2});
figure();
    zplot({Sx,S1,S2},struct('type',3));

if wf > 0
    figure();
        qqplot_(S1,Nf-wf+1);
end
