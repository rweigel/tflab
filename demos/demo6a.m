%% Test of LEMI MT code
% Use LEMI MT program that comes with LEMI hardware to compute TF. Note
% that LEMI MT requires three components of B. (Bz is used to compute
% another TF that is not used here). LEMI MT code is in is available from
% https://github.com/rweigel/lemimt by request. This repository is not 
% public because no license information was found in the LEMI MT software
% package.

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

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
    A.dB = 0.2;
    A.dE = 0.2;
    A.Z  = 1;

Sx_opts = struct('N',N,'n',n,'alpha',alpha,'A',A);
Sx = demo_signals('powerlaw',Sx_opts);

opts = transferfnFD_options(0);
  opts.transferfnFD.loglevel = 0;
  opts.fd.program.name = 'lemimt';
  opts.fd.program.options = ''; % Command line options, e.g., '-r -c'
  opts.info.inunit= 'nT';
  opts.info.outunit= 'mV/km';
  opts.info.timeunit = 's';
  opts.info.timedelta = 1; % time in timeunit between records.
  opts.info.description = 'LEMI';

B(:,1) = Sx.In;
B(:,2) = Sx.In;
B(:,3) = randn(size(B,1),1);

E(:,1) = 0.5*B(:,1) + 0.5*B(:,2);
E(:,2) = 0.5*B(:,1) + 0.5*B(:,2);

S = transferfnFD_lemimt(B,E,opts);
S = transferfnFD_metrics(S,opts);

figure(1);clf;
    tsplot(S,struct('type','raw'));

figure(2);clf;
    tsplot(S,struct('type','error'));

%figure(3);clf;
%    psdplot(S,'raw');

figure(4);clf;
    %zplot(S0a,S1);
    zplot(S);
