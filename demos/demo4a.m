%% Compare OLS and TLS

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

n = 1000; % Number of data points
N = 1000; % Number of frequencies (including 0)

alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 0;
    alpha.Z  = 0;

A = struct();
    A.B  = 1;
    A.dB = 0.2;
    A.dE = 0.2;
    A.Z  = 1;

Sx_opts = struct('N',N,'n',n,'alpha',alpha,'A',A);
Sx = demo_signals(-2,Sx_opts);

opts1 = transferfnFD_options(0);
    opts1.transferfnFD.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[4,4], 'linear'};
    opts1.fd.stack.average.function = '';
    opts1.td.window.width = n/5;
    opts1.td.window.shift = n/5;

S1 = transferfnFD(Sx.In,Sx.Out,opts1);
S1.Options.description = 'Stack Average/OLS';

opts2 = transferfnFD_options(0);
    opts2.transferfnFD.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[4,4], 'linear'};
    opts2.fd.stack.average.function = '';
    opts2.td.window.width = n/5;
    opts2.td.window.shift = n/5;
    opts2.fd.regression.function = @tls_regress;
    opts2.fd.regression.functionstr = 'TLS';
    opts2.fd.regression.functionargs = {};
    opts2.fd.regression.loglevel = 1;

S2 = transferfnFD(Sx.In,Sx.Out,opts2);
S2.Options.description = 'Stack Average/TLS';

% Use same variable labels from S1 for Sa
Sx.Options.info = S1.Options.info; 
Sx.Options.description = 'Actual';

if isfield(Sx,'OutNoise')
    S1.OutNoise = Sx.OutNoise;
    S1.OutNoisePSD = psd(S1.OutNoise);
end
if isfield(Sx,'InNoise')
    S1.InNoise = Sx.InNoise;
    S1.InNoisePSD = psd(S1.InNoise);
end

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
    psdplot({S1,S2},struct('type','error-smoothed'));

figure(f);clf;f = f+1;
    zplot({Sx, S1, S2},struct('type',3));
    
