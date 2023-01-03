%% Single frequency not on DFT grid (so leakage)
% See leakage.m for verifcation that spectrum is not symmetric
% about f when signal phase is not zero or pi/2 (is not pure sin or cos).

clear;

addpath([fileparts(mfilename('fullpath')),'/../plot']);
addpath([fileparts(mfilename('fullpath')),'/..']);


Nt = 100;
f  = 25.5/Nt;
wf = 0;

Sa_opts = struct('Nt',Nt,'Z',1+1j,'f',f,'dB',0.0,'dE',0.0);
Sa = demo_signals(-3,Sa_opts);

opts1 = transferfnFD_options(0);
    opts1.transferfnFD.loglevel = 1;
    opts1.fd.evalfreq.functionargs = {[1,wf], 'linear'};

S1 = transferfnFD(Sa.In,Sa.Out,opts1);
S1.Options.description = 'Estimated';

% Use same variable labels from S1 for Sa
Sa.Options.info = S1.Options.info; 
Sa.Options.description = 'Actual';

set(0,'DefaultFigureWindowStyle','docked')
fn = 1;
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','error'));

figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','error'));

figure(fn);clf;fn = fn+1;
    zplot({Sa,S1},struct('plottype',3));
figure(fn);clf;fn = fn+1;
    zplot({Sa,S1});

if wf > 0    
    figure(fn);clf;fn = fn+1;
        qqplot_(S1,25);
end
