%% Single frequency not on DFT grid (so leakage)
% Demonstrates use of time domain prewhitening (first difference)

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
S1 = tflab(Sx.In,Sx.Out,opts1);
S1.Options.description = 'No prewhiten';

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.prewhiten.function = @prewhiten;
    opts2.td.prewhiten.functionstr = 'First difference';
    opts2.td.prewhiten.functionargs = {'diff'};
S2 = tflab(Sx.In,Sx.Out,opts2);
S2.Options.description = 'Diff prewhiten';

% Use same variable labels from S1 for Sa
Sx.Options.info = S1.Options.info; 
Sx.Options.description = 'Actual';


set(0,'DefaultFigureWindowStyle','docked')
f = 1;
figure(f);clf;f = f+1;
    tsplot(S1,struct('type','raw'));

figure(f);clf;f = f+1;
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure(f);clf;f = f+1;
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure(f);clf;f = f+1;
    psdplot(S1,struct('type','raw'));
figure(f);clf;f = f+1;
    psdplot({S1,S2},struct('type','prewhitened'));    
figure(f);clf;f = f+1;
    psdplot({S1,S2},struct('type','error'));

figure(f);clf;f = f+1;
    zplot({Sx,S1,S2});
figure(f);clf;f = f+1;
    zplot({Sx,S1,S2},struct('type',3));

if wf > 0    
    figure(f);clf;f = f+1;
        qqplot_(S1,Nf-wf+1);
end
