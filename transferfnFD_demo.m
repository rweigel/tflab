clear;

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
    tsplot(S1,struct('type','raw'));   % Plot raw input/output data
figure(2);clf;
    tsplot(S1,struct('type','error')); % Plot raw input/output data
figure(3);clf;
    psdplot(S1,struct('type','raw'));
figure(4);clf;
    zplot({S1,S2});
figure(5);clf;  
    zplot({S1,Sx}); % Compare exact with computed



        
    
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
    tsplot(S1,struct('type','raw'));   % Plot raw input/output data
figure(2);clf;
    tsplot(S1,struct('type','error')); % Plot raw input/output data
%figure(3);clf;
%    psdplot(S1,struct('type','raw'));
figure(4);clf;
    %zplot(S0a,S1);   % Compare exact with computed
    zplot(S1);
%figure(5);clf;
%    hplot(S1,[-5,5]);   % Compare exact with computed



