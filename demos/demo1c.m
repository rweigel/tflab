%% Z nonzero at single frequency on DFT grid (so no leakage)
%  with prewhitening

% The change in the In and Out spectra is due to the fact that the
% pre-whitening filter returns [0; diff(X)]. A zoom-in near t=1 shows
% that harmonic wave has a discontinuity in its slope; the prewhitened
% time series can be decomposed into d*delta(1) + A*sin(2*pi*t/T + phase).
% The d*delta(1) terms results in a flat spectra in the frequency domain.
% To avoid this discontinuity, the pre-whitening filter could return
% a truncated timeseries, diff(X). In this case, both the input and
% output spectrum will be affected by leakage because the truncated 
% timeseries will not have a length that is an integer multiple of T.

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 1000;
k  = 25;   % Frequency index
Z  = (1+1j)/sqrt(2);
f  = k/Nt;
wf = 0;

Sx_opts = struct('Nt', Nt, 'Z',Z, 'f', f, 'dB', 0.0, 'dE', 0.0);
Sx = demo_signals('simple',Sx_opts);
Sx.Options.description = 'Actual';

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
S2.Options.description = 'Diff whitened';


dock on;figure(1);close all;

figure();
    tsplot(S1,struct('type','original'));
figure();
    ax = tsplot(S2,struct('type','whitened'));
    set(ax(1),'XLim',[0,30]);
    set(ax(2),'XLim',[0,30]);

figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure();
    dftplot(S1,struct('type','original-raw'));
figure();
    dftplot(S2,struct('type','original-whitened'));
figure();
    dftplot({S1,S2},struct('type','error'));

figure();
    snplot({S1,S2});

figure();
    zplot({Sx,S1,S2});
figure();
    zplot({Sx,S1,S2},struct('type',3));

if wf > 0 
    figure();
        qqplot_(S1,k);
end
