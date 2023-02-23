%% Flat spectra and frequency dependent phase shift

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
setpaths();

N = 101;
f = fftfreqp(N);
t = (0:N-1)';
for i = 2:length(f)
    A(i,1)   = (i-1);
    Phi(i,1) = -2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t);
    E(:,i) = cos(2*pi*f(i)*t + Phi(i));
    %E(:,i) = cos(2*pi*f(i)*t);
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

figure(1);clf;
    tsplot(S1,struct('type','raw'));   % Plot raw input/output data
figure(2);clf;
    tsplot(S1,struct('type','error')); % Plot raw input/output data
figure(3);clf;
    psdplot(S1,struct('type','raw'));
figure(4);clf;
    zplot(S1);   % Compare exact with computed
figure(5);clf;
    hplot(S1,[-5,5]);   % Compare exact with computed
