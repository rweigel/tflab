%% Single frequency not on DFT grid (so leakage)
% Demonstrates use of time domain prewhitening (first difference)

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
f  = 25.5/Nt;
wf = 0;

Sa_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sa = demo_signals('simple',Sa_opts);
Sa.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
S1 = tflab(Sa.In,Sa.Out,opts1);
S1.Options.description = 'No prewhiten';

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.whiten.function = @whiten;
    opts2.td.whiten.functionstr = 'First difference';
    opts2.td.whiten.functionargs = {'diff'};
S2 = tflab(Sa.In,Sa.Out,opts2);
S2.Options.description = 'Diff prewhiten';

dock on;figure(1);close all;

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S2,struct('type','whitened'));
figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));


figure();
    dftplot(S1,struct('type','original'));
figure();
    dftplot(S2,struct('type','whitened'));
figure();
    dftplot({S1,S2},struct('type','error'));

figure();
    zplot({Sa,S1,S2});
figure();
    zplot({Sa,S1,S2},struct('type',3));

if wf > 0
    figure();
        qqplot_(S1,Nf-wf+1);
end
