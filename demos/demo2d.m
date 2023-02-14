%% Single frequency on DFT grid (so leakage) and zero padding

clear;

addpath([fileparts(mfilename('fullpath')),'/../plot']);
addpath([fileparts(mfilename('fullpath')),'/..']);

Nt = 100;
f  = 25.5/Nt;
wf = 0;
Nz = 2*Nt; % Number of zeros to append.

Sa_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sa = demo_signals(-3,Sa_opts);

opts1 = transferfnFD_options(0);
    opts1.transferfnFD.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};
    % https://dsp.stackexchange.com/a/22041
    % https://dsp.stackexchange.com/a/24426
    opts1.td.zeropad = Nz;
    
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
        qqplot_(S1,Nf-wf+1);
end
