clear;
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=661333
% http://articles.adsabs.harvard.edu//full/1995A%26A...300..707T/0000709.000.html

addpath([fileparts(mfilename('fullpath')),'/plot']);

%% Test effect of adding frequencies in signal between fft grid frequencies
N = 1000; % Number of frequencies
n = 1000; % Number of data points

alpha = struct();
    alpha.B  = [0,   0];
    alpha.dB = [0,   0];
    alpha.dE = 0;
    alpha.Z  = [1, 1];

A = struct();
    A.B  = [1,  1];
    A.dB = [0.0,  0.0];
    A.dE = 0;
    A.Z  = [1,  1];
  
alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 0;
    alpha.Z  = 1;

A = struct();
    A.B  = 1;
    A.dB = 0;
    A.dE = 0;
    A.Z  = 1;

% E(w) = Zo(w)B(w) + dE(w)
% Z(w) = E(w)/B(w)
% Z(w) = Zo(w)B(w)/B(w) + dE(w)/B(w)
% Z(w) = Zo(w) + dE(w)/B(w)

if 0
    B = [];
    E = [];
    for i = 1:10
        %A.dE = 10*i;
        Sx = transferfnFD_demo_signals(-2,struct('N',N,'n',n,'alpha',alpha,'A',A));
        E = [E;Sx.Out];
        B = [B;Sx.In];
    end

    opts = transferfnFD_options(0);
        opts.transferfnFD.loglevel = 1;
        %opts.fd.regression.function = @ols_analytic;    
        %opts.fd.evalfreq.functionargs = {[10,10], 'linear'};
        opts.td.window.width = n;
        opts.td.window.shift = n;
    S1 = transferfnFD(B,E,opts);
end

Sx = transferfnFD_demo_signals(-2,struct('N',10*N,'n',10*n,'alpha',alpha,'A',A));

Estd = std(Sx.Out);
E = reshape(Sx.Out,n,10);
for i = 1:size(E,2)
    E(:,i) = E(:,i) + (Estd)*(i/100)*randn(size(E,1),1);
end
E = E(:);

opts = transferfnFD_options(0);
    opts.transferfnFD.loglevel = 1;
    opts.td.window.width = n;
    opts.td.window.shift = n;
S1 = transferfnFD(Sx.In(1:10*n,:),E,opts);

opts = transferfnFD_options(0);
    opts.transferfnFD.loglevel = 1;
    %opts.fd.evalfreq.functionargs = {[10,10], 'linear'};
    %opts.td.window.width = 10*n;
    %opts.td.window.shift = 10*n;
S2 = transferfnFD(Sx.In,E,opts);

S1.Metrics
S2.Metrics

Sx.Time = S1.Time;                 % Use default computed time from S1 for S0
Sx.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
Sx.Options.description = 'Actual';

close all;
figure(1);clf;
    timeseries_plot(S1,struct('type','raw'));   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,struct('type','error')); % Plot raw input/output data
figure(3);clf;
    spectrum_plot(S1,struct('type','raw'));
figure(4);clf;
    transferfnZ_plot({S1,S2});
figure(5);clf;  
    transferfnZ_plot({S1,Sx}); % Compare exact with computed

break

%% Test effect of adding frequencies in signal between fft grid frequencies
N = 1000; % Number of frequencies
n = 1000; % Number of data points

alpha = struct();
    alpha.B  = [0,   0];
    alpha.dB = [0,   0];
    alpha.dE = 0;
    alpha.Z  = [1, 1];

A = struct();
    A.B  = [1,  1];
    A.dB = [0.0,  0.0];
    A.dE = 0;
    A.Z  = [1,  1];

alpha = struct();
    alpha.B  = 0;
    alpha.dB = 0;
    alpha.dE = 1;
    alpha.Z  = 1;

A = struct();
    A.B  = 1;
    A.dB = 0;
    A.dE = 0.1;
    A.Z  = 1;
    
Sx = transferfnFD_demo_signals(-2,struct('N',10*N,'n',10*n,'alpha',alpha,'A',A));

opts = transferfnFD_options(1);
    opts.transferfnFD.loglevel = 1;
    %opts.fd.regression.function = @ols_analytic;    
    %opts.fd.evalfreq.functionargs = {[1,1], 'linear'};
    opts.td.window.width = n;
    opts.td.window.shift = n;
S1 = transferfnFD(Sx.In,Sx.Out,opts);

opts = transferfnFD_options(1);
    opts.transferfnFD.loglevel = 1;
    %opts.fd.evalfreq.functionargs = {[1,1], 'linear'};
    opts.td.window.width = 10*n;
    opts.td.window.shift = 10*n;
S2 = transferfnFD(Sx.In,Sx.Out,opts);

if 0
    B(:,3) = randn(size(B,1),1);
    E(:,2) = E(:,1);
    S1 = transferfnFD_lemimt(B,E,opts);

    S1.In = B(:,1:2);
    S1.Z = S1.Z(:,1:2);
    S1.Out = E(:,1:1);
    S1.Time = [0:size(B,1)-1]';
    S1.Options = opts;
    S1 = transferfnMetrics(S1, opts);
end

Sx.Time = S1.Time;                 % Use default computed time from S1 for S0
Sx.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
Sx.Options.description = 'Actual';

close all;
figure(1);clf;
    timeseries_plot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,'error'); % Plot raw input/output data
figure(3);clf;
    spectrum_plot(S1,'raw');
figure(4);clf;
    transferfnZ_plot(S1,S2);   
figure(5);clf;  
    transferfnZ_plot(S1,Sx); % Compare exact with computed

break    
    loglog(1./S0.fe(2:end),abs(S0.Z(2:end,1,:)),'k','LineWidth',5);
    hold on;grid on;
    loglog(1./S2.fe(2:end),abs(S2.Z(2:end,1,:)),'b','LineWidth',3);
    loglog(1./S1.fe(2:end),abs(S1.Z(2:end,1,:)),'r','LineWidth',3);
    %loglog(1./S1.fe(2:end),abs(squeeze(S1.Segment.Z(2:end,1,:))))
    
        
break
    
%% Test of lemimt code
N = 10001; % Number of frequencies
n = 1001;  % Number of data points

S0a = transferfnFD_demo_signals(-2,struct('N',N,'n',n));
S0b = transferfnFD_demo_signals(-2,struct('N',N,'n',n));

opts = transferfnFD_options(0);
addpath('~/git/lemimt');
%opts.fd.program.name = 'lemimt'
opts.transferfnFD.loglevel = 0;

%B(:,1) = randn(N,1);
%B(:,2) = randn(N,1);
B(:,1) = S0a.In;
B(:,2) = S0b.In;
B(:,3) = randn(size(B,1),1);

E(:,1) = 0.5*B(:,1) + 0.5*B(:,2);
E(:,2) = 0.5*B(:,1) + 0.5*B(:,2);
%E(:,3) = zeros(N,1);
%E(:,4) = zeros(N,1);

if 0
S1 = transferfnFD_lemimt(B,E,opts);

S1.In = B(:,1:2);
S1.Z = S1.Z(:,1:2);
S1.Out = E(:,1:1);
S1.Time = [0:size(B,1)-1]';
S1.Options = opts;
S1 = transferfnMetrics(S1, opts);
end

S1 = transferfnFD(B,E,opts);

figure(1);clf;
    timeseries_plot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    timeseries_plot(S1,'error'); % Plot raw input/output data
%figure(3);clf;
%    spectrum_plot(S1,'raw');
figure(4);clf;
    %transferfnZ_plot(S0a,S1);   % Compare exact with computed
    transferfnZ_plot(S1);
%figure(5);clf;
%    transferfnH_plot(S1,[-5,5]);   % Compare exact with computed


break


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
