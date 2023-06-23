%% Test of LEMI MT code
% Use LEMI MT program that comes with LEMI hardware to compute TF. Note
% that LEMI MT requires three components of B. (Bz is used to compute
% another TF that is not used here). LEMI MT code is in is available from
% https://github.com/rweigel/lemimt by request. This repository is not 
% public because no license information was found in the LEMI MT software
% package.

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

lemimt_dir = '~/git/lemimt';
if ~exist(lemimt_dir,'dir')
    url = 'https://github.com/rweigel/lemimt';
    error('Need to download and compile lemimt code from %s', url);
end
addpath(lemimt_dir);

n = 1000; % Number of data points
N = 1000; % Number of frequencies (including 0)

alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 0;
    alpha.Z  = 0;

% Amplitudes
A = struct();
    A.B  = 1;
    A.E  = 1;
    A.dB = 0.1;
    A.dE = 0.1;
    A.Z  = 1;

% Here we set rng(1) so that same In and Out are created each time this
% script is run. This is to avoid error that sometimes appears:
% "Note: The following floating-point exceptions are signalling:
% IEEE_INVALID_FLAG"
% TODO: Determine what is causing this error.
rng(1);
Sx_opts = struct('N',N,'n',n,'alpha',alpha,'A',A);
Sx = demo_signals('powerlaw',Sx_opts);

opts = tflab_options(1);
  opts.tflab.loglevel = 0;
  opts.fd.program.name = 'lemimt';
  opts.fd.program.options = ''; % Command line options, e.g., '-r -c'

B(:,1) = Sx.In;
E(:,1) = Sx.Out;

S = tflab_lemimt(B,E);

S.Options.description = 'LEMI MT';
S.Metadata.inunit= 'nT';
S.Metadata.outunit= 'mV/km';
S.Metadata.timeunit = 's';
S.Metadata.timedelta = 1; % time in timeunit between records.

% Set options used to computed averaged coherence and signal to error.
S.Options.fd.evalfreq.function = @evalfreq;
S.Options.fd.evalfreq.functionargs = {7, 'logarithmic'};
S.Options.fd.window.function = @rectwin;

S = tflab_metrics(S);

dockreset();

figure();clf;
    tsplot(S,struct('type','original'));

figure();clf;
    tsplot(S,struct('type','error'));

figure();clf;
    %zplot(S0a,S1);
    zplot(S);
