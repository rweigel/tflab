%% Z nonzero at single frequency on DFT grid (so no leakage) and zero padding
% See also
%   Motivation for zero padding: https://dsp.stackexchange.com/a/22041
%   Proof that zero padding corresponds to sinc-weighted average in
%   frequency domain: https://dsp.stackexchange.com/a/24426 (in the post
%   they use "interpolation" instead of "weighted average").

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

Nt = 100;
k  = 25;
Z  = (1+1j)/sqrt(2);
wf = 1;

regstr = sprintf('OLS/$N_b=%d$',2*wf+1); 

Sa_opts = struct('Nt',Nt,'Z',Z,'k',k,'dB',0,'dE',0);
Sa = demo_signals('simple',Sa_opts);
Sa.Options.description = 'Actual';

opts1 = tflab_options(0);
    opts1.tflab.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.zeropad = Nt;
    
S1 = tflab(Sa.In,Sa.Out,opts1);
S1.Options.description = sprintf('zeropadded/%s',regstr);

dockreset();

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','error'));
    
figure();
    dftplot(S1,struct('type','original'));
figure();
    dftplot(S1,struct('type','zeropadded'));    
figure();
    dftplot(S1,struct('type','error'));

figure();
    zplot({Sa,S1});
figure();
    zplot({Sa,S1},struct('type',3));

if wf > 0    
    figure();
        [~,Ik] = min(abs(1./S1.fe-Nt/k));
        qqplot_(S1,struct(),Ik);
end