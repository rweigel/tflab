%% Average of 5 vs full time series

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

N = 1000; % Number of frequencies
n = 1000; % Number of data points

alpha = struct();
    alpha.B  = [0,   0];
    alpha.dB = [0,   0];
    alpha.dE = 0;
    alpha.Z  = [1, 1];

A = struct();
    A.B  = [1,  1];
    A.dB = [0.0,  0.0];
    A.dE = 0;
    A.Z  = [1,  1];

alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 1;
    alpha.Z  = 1;

A = struct();
    A.B  = 1;
    A.dB = 0;
    A.dE = 0.1;
    A.Z  = 1;
    
Sx = demo_signals('powerlaw',struct('N',10*N,'n',10*n,'alpha',alpha,'A',A));

opts = tflab_options(1);
    opts.tflab.loglevel = 1;
    %opts.fd.regression.function = @ols_analytic;    
    %opts.fd.evalfreq.functionargs = {[1,1], 'linear'};
    opts.td.window.width = n;
    opts.td.window.shift = n;
S1 = tflab(Sx.In,Sx.Out,opts);

opts = tflab_options(1);
    opts.tflab.loglevel = 1;
    %opts.fd.evalfreq.functionargs = {[1,1], 'linear'};
    opts.td.window.width = 10*n;
    opts.td.window.shift = 10*n;
S2 = tflab(Sx.In,Sx.Out,opts);

if 0
    B(:,3) = randn(size(B,1),1);
    E(:,2) = E(:,1);
    S1 = tflab_lemimt(B,E,opts);

    S1.In = B(:,1:2);
    S1.Z = S1.Z(:,1:2);
    S1.Out = E(:,1:1);
    S1.Time = [0:size(B,1)-1]';
    S1.Options = opts;
    S1 = tflab_metrics(S1, opts);
end

Sx.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
Sx.Options.description = 'Actual';

dock on;figure(1);close all;

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','error'));
figure();
    dftplot(S1,struct('type','original'));
figure();
    zplot({S1,S2});   
figure();
    zplot({S1,Sx}); % Compare exact with computed
