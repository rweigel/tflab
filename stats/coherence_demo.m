clear
%% Demo 1
N = 256*4;
f = fftfreqp(N/4); % Get positive DFT frequencies for a signal of length N
t = (0:N-1)';    % Time index

% Create input/output signals.

E = zeros(N,1);
B = zeros(N,1);
i = 32;
% Uniformly distributed random phase in [-pi,pi]
p = 0*pi*(-1 + 2*rand(1)); 
B = B + cos(2*pi*f(i)*t + p);             % Input
%E = E + cos(2*pi*f(i)*t + p + 2*pi*f(i)); % Output
E = E + cos(2*pi*f(i)*t + p + pi/2); % Output
B = B + 0.2*randn(N,1); % Add input noise
E = E + 0.2*randn(N,1); % Add output noise

[Cxy,fxy] = mscohere(E,B,rectwin(N/4),0,N/4);

opts = tflab_options(0);
    opts.tflab.loglevel = 1;
    opts.fd.evalfreq.functionargs = {[1,5], 'linear'};
[Cxy2,fxy2] = coherence(E,B,1,opts);

% Manual version of mscohere with given arguments
[dftE,fxy3] = fftu(reshape(E,256,4));
[dftB,~] = fftu(reshape(B,256,4));
sxx = abs(sum(dftE.*conj(dftE),2));
syy = abs(sum(dftB.*conj(dftB),2));
sxy = abs(sum(dftE.*conj(dftB),2));
Cxy3 = (sxy.^2)./(sxx.*syy);

figure(1);clf
    plot(E);
    hold on;
    plot(B);
figure(2);clf;
    plot(fxy/(2*pi),Cxy,'b','linewidth',4);
    hold on;
    plot(fxy2,Cxy2,'r','LineWidth',2)
    plot(fxy3,Cxy3,'g','LineWidth',2)
    legend('mscohere','coherence','manual welch');

