%% Single frequency not on DFT grid (so leakage)
% See notes/leakage.m for additional demos of leakage and comparison
% with analytical formulas. In particular, it shows that spectrum is not
% symmetric about f when signal phase is not zero or pi/2 (is not pure
% sin or cos).

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
f  = 25.5/Nt;
wf = 0;
dB = 0.0;
dE = 0.0;
Z  = (1+1j)/sqrt(2);

regstr = sprintf('OLS/$N_b=%d$',2*wf+1); 

Sa_opts = struct('Nt',Nt,'Z',Z,'f',f,'dB',dB,'dE',dE);
Sa = demo_signals('simple',Sa_opts);
Sa.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};

S1 = tflab(Sa.In,Sa.Out,opts1);
S1.Options.description = sprintf('%s/dIn=%.1f/dOut=%0.1f',...
                                 regstr, dB, dE);

dockreset();

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','error'));

figure();
    dftplot(S1,struct('type','original'));
figure();
    dftplot(S1,struct('type','error'));

figure();
    snplot(S1);

figure();
    zplot({Sa,S1},struct('type',3));
figure();
    zplot({Sa,S1});

if wf > 0    
    figure();
        qqplot_(S1,25);
end
