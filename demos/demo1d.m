%% Z nonzero at single frequency on DFT grid (so no leakage) and zero padding
% See also
%   Motivation for zero padding: https://dsp.stackexchange.com/a/22041
%   Proof that zero padding corresponds to sinc-weighted average in
%   frequency domain: https://dsp.stackexchange.com/a/24426 (in the post
%   they use "interpolation" instead of "weighted average").

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

Nt = 100;
k  = 25;
wf = 0;

Sa_opts = struct('Nt',Nt,'Z',1j,'k',k,'dB',0,'dE',0);
Sa = demo_signals('simple',Sa_opts);

opts1 = transferfnFD_options(0);
    opts1.transferfnFD.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    opts1.td.zeropad = Nt;
    
S1 = transferfnFD(Sa.In,Sa.Out,opts1);

S1.Options.description = 'Estimated';

% Use same variable labels from S1 for Sa
Sa.Options.info = S1.Options.info; 
Sa.Options.description = 'Actual';

set(0,'DefaultFigureWindowStyle','docked')
f = 1;
figure(f);clf;f = f+1;
    tsplot(S1,struct('type','raw'));
figure(f);clf;f = f+1;
    tsplot(S1,struct('type','error'));
    
figure(f);clf;f = f+1;
    psdplot(S1,struct('type','raw'));
figure(f);clf;f = f+1;
    psdplot(S1,struct('type','zeropadded'));    
figure(f);clf;f = f+1;
    psdplot(S1,struct('type','error'));

figure(f);clf;f = f+1;
    zplot({Sa,S1});
figure(f);clf;f = f+1;
    zplot({Sa,S1},struct('type',3));

if wf > 0    
    figure(f);clf;f = f+1;
        qqplot_(S1,k);
end