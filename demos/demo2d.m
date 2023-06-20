%% Single frequency on DFT grid (so leakage) and zero padding

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
f  = 25.0/Nt;
wf = 1;
Nz = Nt; % Number of zeros to append.

Sa_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sa = demo_signals('simple',Sa_opts);
Sa.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
S1 = tflab(Sa.In,Sa.Out,opts1);
S1.Options.description = 'No padding';

opts2 = tflab_options(0);
    opts2.tflab.loglevel = 1;
    opts2.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts2.td.zeropad = Nz;

S2 = tflab(Sa.In,Sa.Out,opts2);
S2.Options.description = sprintf('Padded with %d zeros',Nz);

dock on;figure(1);pause(0.1);close all;

figure();
    tsplot(S1,struct('type','original'));

figure();
    tsplot(S2,struct('type','zeropadded'));
    
figure();
    tsplot(S1,struct('type','error','title',S1.Options.description));
figure();
    tsplot(S2,struct('type','error','title',S2.Options.description));

figure();
    dftplot(S2,struct('type','original-raw-reals'));    
figure();
    dftplot(S2,struct('type','original-raw-imaginaries'));    
figure();
    dftplot(S2,struct('type','zeropadded-raw-reals'));    
figure();
    dftplot(S2,struct('type','zeropadded-raw-imaginaries'));    
figure();
    dftplot({S1,S2},struct('type','error'));

figure();
    zplot({Sa,S1,S2});
figure();
    zplot({Sa,S1,S2},struct('type',3));

if wf > 0
    figure();
        [~,fidx] = min(abs(S1.fe-f));
        qqplot_(S1,fidx);
    figure();
        [~,fidx] = min(abs(S2.fe-f));
        qqplot_(S2,fidx);        
end