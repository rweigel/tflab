%% Z nonzero at single frequency on DFT grid (so no leakage)

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;  % Number of time samples
k  = 10;   % Freq. index of non-zero Z
Z  = (1+1j)/sqrt(2);
wf = 0;    % 2*wf + 1 is number of points for regression.

Sa_opts = struct('Nt',Nt,'Z',Z,'k',k,'dB',0,'dE',0);
tfa = demo_signals('simple',Sa_opts);
tfa.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};

tf1 = tflab(tfa.In, tfa.Out, opts1);

dockreset();

figure();
    tsplot(tf1,struct('type','original'));
figure();
    tsplot(tf1,struct('type','error'));

figure();
    dftplot(tf1,struct('type','original-raw-magnitudes'));
figure();
    dftplot(tf1,struct('type','original-raw-phases'));
figure();
    dftplot(tf1,struct('type','error-raw-magphase'));

figure();
    snplot(tf1);
    
figure();
    zplot({tfa,tf1});
figure();
    zplot({tfa,tf1},struct('type',3));

if wf > 0
    % Regression was used. Show qq plot.
    figure();
        qqplot_(tf1,k);
end
