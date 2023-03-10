clear

addpath(fullfile(fileparts(mfilename('fullpath'))));
tflab_setpaths();

%% Demo 1
% Input:  1-D time series with a flat spectrum and random phase.
% Output: 1-D time series with a linear frequency dependent phase.

N = 99;
f = fftfreqp(N); % Get positive DFT frequencies for a signal of length N
t = (0:N-1)';    % Time index

% Create input/output signals.

E = zeros(N,1);
B = zeros(N,1);
for i = 2:length(f)
    % Uniformly distributed random phase in [-pi,pi]
    p = pi*(-1 + 2*rand(1)); 
    B = B + cos(2*pi*f(i)*t + p);             % Input
    E = E + cos(2*pi*f(i)*t + p + 2*pi*f(i)); % Output
    B = B + 0.0*randn(N,1); % Add input noise
    E = E + 0.0*randn(N,1); % Add output noise
end

% Remove mean
B = B - mean(B);
E = E - mean(E);

% Scale to [-1,1]
B = B/max(B);
E = E/max(E);

% Set means
E = E + 2.0;
B = B + 1.0;

% Get default algorithm options
opts = tflab_options(0);
printstruct(opts);

% Estimate Z
tf1 = tflab(B,E,opts);

vs_period = 0; % (Default is to plot vs. period.)

set(0,'DefaultFigureWindowStyle','docked')

fn = 1;
figure(fn);clf;fn = fn+1;
    tsplot(tf1,struct('type','raw','vs_period',vs_period));
figure(fn);clf;fn = fn+1;
    tsplot(tf1,struct('type','error','vs_period',vs_period));

figure(fn);clf;fn = fn+1;
    psdplot(tf1,struct('type','raw','vs_period',vs_period));
figure(fn);clf;fn = fn+1;
    psdplot(tf1,struct('type','raw-phase','vs_period',vs_period));
figure(fn);clf;fn = fn+1;
    psdplot(tf1,struct('type','error','vs_period',vs_period));

figure(fn);clf;fn = fn+1;
    zplot(tf1,struct('vs_period',vs_period));
figure(fn);clf;fn = fn+1;
    zplot(tf1,struct('type',3,'vs_period',vs_period));
