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
    opts1.description = 'Estimated';

tf1 = tflab(tfa.In, tfa.Out, opts1);

dock on;figure(1);close all;

figure();
    tsplot(tf1,struct('type','raw'));
figure();
    tsplot(tf1,struct('type','error'));

figure();
    psdplot(tf1,struct('type','raw'));
figure();
    psdplot(tf1,struct('type','raw-phase'));
figure();
    psdplot(tf1,struct('type','error'));

figure();
    snplot(tf1);
    
figure();
    zplot({tfa,tf1});
figure();
    zplot({tfa,tf1},struct('type',3));

if wf > 0
    % Regression was used. Show qq plot.
    figure(fn);clf;fn = fn+1;
        qqplot_(tf1,k);
end


%% Compute variance of Z as a function of Nt

if 0
    Ikeep = 1:size(tf1.Z,1);
    Ikeep = Ikeep(Ikeep ~= k);
    Z = tf1.Z(Ikeep);
    fprintf('Nt = %f; var(Z) = %f\n\n',Nt,var(Z(2:end-1)));
    
    opts1 = tflab_options(0);
        opts1.tflab.loglevel = 0;
        opts1.fd.evalfreq.functionargs = {[1,0], 'linear'};

    for Nt = [10^2,10^3,10^4,10^5] % [2^10,2^12,2^14]

        k = round(N/4);

        Sa_opts = struct('Nt',Nt,'Z',1+1j,'k',k,'dB',0,'dE',0);
        tfa = demo_signals('simple',Sa_opts);

        tf1 = tflab(tfa.In,tfa.Out,opts1);
        Ikeep = 1:size(tf1.Z,1);
        Ikeep = Ikeep(Ikeep ~= k+1);
        Z = tf1.Z(Ikeep);
        fprintf('Nt = %f; var(Z) = %f\n',Nt,var(Z(2:end-1)));
    end
    % Nt = 100.000000; var(Z) = 2.305222
    % Nt = 1000.000000; var(Z) = 30.523350
    % Nt = 10000.000000; var(Z) = 15.500834
    % Nt = 100000.000000; var(Z) = 10.954189
end
