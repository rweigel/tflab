clear

addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

%% Demo 1
% Input signal a 1-D time series that has flat spectrum and random phase.
% (If phase is not random, signal is delta function.)
% Output signal is a 1-D time series with a linear frequency dependent phase.

N = 101;
f = fftfreqp(N); % Get positive DFT frequencies for a signal of length N
t = (0:N-1)';    % Time index

% Create input/output signals.

E = zeros(N,1);
B = zeros(N,1);
for i = 2:length(f)
    phio = pi*(-1 + 2*rand(1));         % Random phase in [-pi,pi]
    B = B + cos(2*pi*f(i)*t + phio);    % Input
    dE = 0.0*randn(N,1);
    E  = dE + E + cos(2*pi*f(i)*t + phio - 2*pi*f(i)); % Output
end

% Scale to [-1,1]
B = B/max(B);
E = E/max(E);

% Get default options
opts = transferfnFD_options(0); 

% Compute estimate of Z
S1 = transferfnFD(B,E,opts);

set(0,'DefaultFigureWindowStyle','docked')
fn = 1;
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    tsplot(S1,struct('type','error'));

figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','raw'));
figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','raw-phase'));
figure(fn);clf;fn = fn+1;
    psdplot(S1,struct('type','error'));

figure(fn);clf;fn = fn+1;
    zplot(S1);
figure(fn);clf;fn = fn+1;
    zplot(S1,struct('type',3));
