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
    error('Need to download and compile lemimt code from https://github.com/rweigel/lemimt');
end
addpath(lemimt_dir);

N = 10001; % Number of frequencies
n = 1001;  % Number of data points

S0a = demo_signals(-2,struct('N',N,'n',n));
S0b = demo_signals(-2,struct('N',N,'n',n));

opts = transferfnFD_options(0);
  opts.fd.program.name = 'lemimt';
  opts.transferfnFD.loglevel = 0;

B(:,1) = S0a.In;
B(:,2) = S0b.In;
B(:,3) = randn(size(B,1),1);

E(:,1) = 0.5*B(:,1) + 0.5*B(:,2);
E(:,2) = 0.5*B(:,1) + 0.5*B(:,2);

S = transferfnFD_lemimt(B,E);
S.In = B(:,1:2);
S.Out = E;
S = transferfnFDMetrics(S{tf},opts2,S{1}.Segment.IndexRange);
S.Options.description = 'LEMI';

figure(1);clf;
    tsplot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    tsplot(S1,'error'); % Plot raw input/output data
%figure(3);clf;
%    psdplot(S1,'raw');
figure(4);clf;
    %zplot(S0a,S1);   % Compare exact with computed
    zplot(S1);
%figure(5);clf;
%    hplot(S1,[-5,5]);   % Compare exact with computed
