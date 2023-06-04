%% Flat spectra and frequency dependent phase shift

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

N = 11;
f = fftfreqp(N);
t = (0:N-1)';
for i = 2:length(f)
    A(i,1)   = (i-1);
    Phi(i,1) = -2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t);
    E(:,i) = cos(2*pi*f(i)*t + Phi(i));
end

B = sum(B,2);
E = sum(E,2);

B = B/max(B);
E = E/max(E);

% Force E to have same mean as B so that H is not offset.
E = E + (mean(B)-mean(E));

% Compute estimate of Z
opts = tflab_options(0);
S1 = tflab(B,E,opts);

dock on;figure(1);close all;

figure();
    tsplot(S1,struct('type','original'));
figure();
    tsplot(S1,struct('type','error'));
figure();
    dftplot(S1,struct('type','original'));
figure();
    zplot(S1);
figure();
    hplot(S1,[-5,5]);   % Compare exact with computed
