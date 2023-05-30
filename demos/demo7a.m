addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;  % Number of time samples
k  = 10;   % Freq. index of non-zero Z
Z  = (1+1j)/sqrt(2);
wf = 0;    % 2*wf + 1 is number of points for regression.

Sa_opts = struct('Nt',Nt,'Z',Z,'k',k,'dB',0,'dE',0);
tfa = demo_signals('simple',Sa_opts);
tfa.Options.description = 'Actual';
n = size(tfa.In,1)-1;
tfa.In = tfa.In + (0:n)'/n;

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.description = 'Estimated';
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.detrend.function = @bandpass_;
    opts1.td.detrend.functionstr = 'Bandpass';  % Optional descriptive name
    opts1.td.detrend.functionargs = {[1/20,0.5]}; % Arguments after first argument to fn.

    
tf1 = tflab(tfa.In, tfa.Out, opts1);


dock on;figure(1);close all;

figure();clf;
    ax = tsplot(tf1,struct('type','raw'));
figure();clf;
    tsplot(tf1,struct('type','Detrended','title','detrended'));
