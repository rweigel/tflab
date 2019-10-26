clear;
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=661333
% http://articles.adsabs.harvard.edu//full/1995A%26A...300..707T/0000709.000.html

addpath([fileparts(mfilename('fullpath')),'/plot']);

%% Test of lemimt code
N = 10001; % Number of frequencies
n = 1001;  % Number of data points

S0 = transferfnFD_demo_signals(-2,struct('N',N,'n',n));

opts = transferfnFD_options(0);
addpath('~/git/lemimt');
opts.fd.program.name = 'lemimt'
opts.transferfnFD.loglevel = 0;

S1 = transferfnFD(S0.In,S0.Out,opts);

break

%% Test effect of adding frequencies in signal between fft grid frequencies
N = 10001; % Number of frequencies
n = 1001;  % Number of data points

S0 = transferfnFD_demo_signals(-2,struct('N',N,'n',n));

opts = transferfnFD_options(0);
opts.transferfnFD.loglevel = 0;

S1 = transferfnFD(S0.In,S0.Out,opts);

S0.Time = S1.Time;                 % Use default computed time from S1 for S0
S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
S0.Options.description = 'Actual';

figure(1);clf;
    timeseries_plot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,'error'); % Plot raw input/output data
figure(3);clf;
    spectrum_plot(S1,'raw');
figure(4);clf;
    transferfnZ_plot(S0,S1);   % Compare exact with computed
    %transferfnZ_plot(S1);
figure(5);clf;
    transferfnH_plot(S1,[-5,5]);   % Compare exact with computed

break

%% Test of E(t+1) = B(t)
N = 11;
H = [0,1]';

S0 = transferfnFD_demo_signals(-1,struct('H',H,'N',N));

opts = transferfnFD_options(0);
opts.transferfnFD.loglevel = 0;
S1 = transferfnFD(S0.In,S0.Out,opts);

S0.Time = S1.Time;                 % Use default computed time from S1 for S0
S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
S0.Options.description = 'Actual';

figure(1);clf;
    timeseries_plot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,'error'); % Plot raw input/output data
figure(3);clf;
    spectrum_plot(S1,'raw');
figure(4);clf;
    transferfnZ_plot(S0,S1);   % Compare exact with computed
figure(5);clf;
    transferfnH_plot(S1,[-5,5]);   % Compare exact with computed

assert(all(abs(real(S1.Z) - real(S0.Z))) <= eps);
assert(all(abs(imag(S1.Z) - imag(S0.Z))) <= eps);

break


clear E B
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
opts = transferfnFD_options(0);
S1 = transferfnFD(B,E,opts);

figure(1);clf;
    timeseries_plot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,'error'); % Plot raw input/output data
figure(3);clf;
    spectrum_plot(S1,'raw');
figure(4);clf;
    transferfnZ_plot(S1);   % Compare exact with computed
figure(5);clf;
    transferfnH_plot(S1,[-5,5]);   % Compare exact with computed



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output signal using given H
H = [1,0,0]';
N = 1000; % Length of input/output signals
S0 = transferfnFD_demo_signals(0,struct('H',H,'N',N));

% Compute estimate of Z
opts = transferfnFD_options(0);
S1 = transferfnFD(S0.In,S0.Out,opts);

S0.Time = S1.Time;                 % Use default computed time from S1 for S0
S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
S0.description = 'Actual';

timeseries_plot(S1,'raw'); % Plot raw input/output data
transferfnZ_plot(S0,S1);   % Compare exact with computed
transferfnH_plot(S0,S1,[-10,10]);   % Compare exact with computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
