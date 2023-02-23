%% Comparison of Stack Regression with Stack Averaging

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

n = 1000; % Number of data points
N = 1000; % Number of frequencies (including 0)
wf = 1;
W  = n/5;

alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 0;
    alpha.Z  = 0;

A = struct();
    A.B  = 1;
    A.dB = 0.5;
    A.dE = 0.0;
    A.Z  = 1;

Sx_opts = struct('N',N,'n',n,'alpha',alpha,'A',A);
Sx = demo_signals('powerlaw',Sx_opts);

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.fd.stack.average.function = '';
    opts1.td.window.width = W;
    opts1.td.window.shift = W;

S1 = tflab(Sx.In,Sx.Out,opts1);
S1.Options.description = 'Stack Regression';

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.window.width = W;
    opts2.td.window.shift = W;

S2 = tflab(Sx.In,Sx.Out,opts2);
S2.Options.description = 'Stack Average';

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
fn = 1;
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure(fn);clf;fn = fn+1;
    tsplot(S2,struct('type','error','title',S2.Options.description));
    
figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    psdplot({S1,S2},struct('type','error-smoothed'));

figure(fn);clf;fn = fn+1;
    zplot({Sx, S1, S2},struct('type',3));
    
if wf > 0
    figure(fn);clf;fn = fn+1;
        qqplot_(S1,2);
end
