clear;
addpath('../../time/')
addpath('../../stats/')
addpath('../misc/')

writeimgs = 0;

tau  = 10;  % Filter decay constant
Ntau = 10;  % Number of filter coefficients = Ntau*tau + 1
N    = 1e4; % Simulation length
Nc   = 50;  
Nf   = 2;
Nss  = 4;   % Ntau*tau*Nss

% IRF for dx/dt + x/tau = delta(0), and ICs
% x_0 = 0 dx_0/dt = 0 approximated using forward Euler.
dt = 1;
gamma = (1-dt/tau);
for i = 1:Ntau*tau
    h1(i,1) = gamma^(i-1);
    h2(i,1) = gamma^(i-1);
    t(i,1) = dt*(i-1);
end
t   = [0:length(h1)]';

%h1  = exp(-(t(2:end)/3/tau).^2);
%h2  = h1;

h1  = [0;h1];
h2  = [0;h2];
hxx =  h1;
hxy = -h1;
hyx =  h2/2;
hyy = -h2/2;

% Create signals
B(:,1) = randn(N,1);
B(:,2) = randn(N,1);

E(:,1) = filter(hxx,1,B(:,1)) + filter(hxy,1,B(:,2));
E(:,2) = filter(hyx,1,B(:,1)) + filter(hyy,1,B(:,2));

E = E(1:end,:);
B = B(1:end,:);

if (1)
for i = 1:size(B,2)
    B(:,i)  = B(:,i) - mean(B(:,i));
end
for i = 1:size(E,2)
    E(:,i) = E(:,i) - mean(E(:,i));
end
end

% Remove non-steady-state part of signals
B  = B(Nss*Nc+1:end,:);
E  = E(Nss*Nc+1:end,:);

fh = [0:Nc/2]'/Nc;

[ZR,feR,HR,tR] = transferfnTD2(B,E,Nc,0);

if (0)
    [ZR,feR] = transferfnFD(B,E,2,'rectangular',Nf);
    ZRi      = Zinterp(feR,ZR,fh);
    HR       = Z2H(feR,ZR);
    HR       = fftshift(HR,1);
    NR       = (size(HR,1)-1)/2;
    tR       = [-NR:NR]';
    EpRH     = Hpredict(tR,HR,B);
    EpRZ     = Zpredict(feR,ZR,B);
end

set(0,'DefaultFigureWindowStyle','docked');
figure(1);clf;hold on;grid on;
    plot(t,hxx,'k.','MarkerSize',20);
    plot(t,hxy,'k.','MarkerSize',20);
    plot(t,hyx,'k.','MarkerSize',20);
    plot(t,hyy,'k.','MarkerSize',20);
    plot(tR,HR(:,1),'r.','MarkerSize',10);
    plot(tR,HR(:,2),'g.','MarkerSize',10);
    plot(tR,HR(:,3),'b.','MarkerSize',10);
    plot(tR,HR(:,4),'c.','MarkerSize',10);
