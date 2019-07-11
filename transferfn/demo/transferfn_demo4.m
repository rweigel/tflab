clear

addpath('../../time/')
addpath('../../stats/')
addpath('../misc/')

writeimgs = 1;

tau = 10;  % Filter decay constant
N   = 86400*10; % Simulation length
Nc  = tau*20;  % Comment out to use Nc = length(h)
df  = 50;  % Width of rectangualar window
nb  = 0.0; % Noise in B
ne  = 0.1; % Noise in E
ndb = 0.0; % Noise in dB

paramstring = sprintf('_ne_%.1f',ne);

% IRF for dx/dt + x/tau = delta(0), and ICs
% x_0 = 0 dx_0/dt = 0 approximated using forward Euler.
dt = 1;
gamma = (1-dt/tau);
for i = 1:10*tau
    h(i,1) = gamma^(i-1);
    t(i,1) = dt*(i-1);
end
hstr = sprintf('(1-1/%d)^{t}; t=1 ... %d; h_{xy}(0)=0;', tau, length(h));
% Add a zero because MATLAB filter function requires it.  Also,
% not having it causes non-physical phase drift with frequency.
h  = [0;h];
t  = [0:length(h)-1]';

% Exact transfer function
Zxy = fft(h);
Nh  = length(Zxy);
Zxy = Zxy(1:floor(Nh/2)+1);
fh  = [0:floor(Nh/2)]'/Nh;

% Exact transfer Function Phase
Pxy = (180/pi)*atan2(imag(Zxy),real(Zxy));

% Add extra values to get nice length
% (because we cut off non-steady state).
N = N + 2*length(h);

% Noise
NE  = [ne*randn(N,1),ne*randn(N,1)];
NB  = [nb*randn(N,1),nb*randn(N,1)];
NdB = [ndb*randn(N,1),ndb*randn(N,1)];

% Create signals
B(:,1) = randn(N,1);
B(:,2) = randn(N,1);

E(:,2) = filter(h,1,B(:,1)+NB(:,1)) + NE(:,2);
E(:,1) = filter(h,1,B(:,2)+NB(:,2)) + NE(:,1);

dB = diff(B);
dB = [dB;dB(end,:)];

% Remove non-steady-state part of signals
B  = B(2*length(h)+1:end,:);
E  = E(2*length(h)+1:end,:);
dB = dB(2*length(h)+1:end,:);

NE  = NE(2*length(h)+1:end,:);
NB  = NB(2*length(h)+1:end,:);
NdB = NdB(2*length(h)+1:end,:);

N = size(B,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time domain
tic()
[Z1,fe1,H1,t1,Ep1] = transferfnTD(B,E,Nc);
fprintf('Time for transferfnTD:    %.2f\n',toc())

tic();
Na = 2;
[Z2,f2,H2,t2,Ep2] = transferfnTDri(B,E,Nc,0,Na);
fprintf('Time for transferfnTDri:  %.2f\n',toc())

tic();
Ns = 86400/10;
[Z3,f3,H3,t3,Ep3] = transferfnTDave(B,E,Nc,0,Ns);
fprintf('Time for transferfnTDave: %.2f\n',toc())

clf;hold on;grid on;grid minor;
    plot(t,h,'k.','MarkerSize',30);
    plot(t1,H1(:,2),'g.','MarkerSize',20);
    plot(t2,H2(:,2),'m.','MarkerSize',15);
    plot(t3,H3(:,2),'b.','MarkerSize',10);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
