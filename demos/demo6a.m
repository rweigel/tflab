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
S1 = transferfnFD_metrics(S1, opts);
end

S1 = transferfnFD(B,E,opts);

figure(1);clf;
    tsplot(S1,'raw');   % Plot raw input/output data
figure(2);clf;
    tsplot(S1,'error'); % Plot raw input/output data
%figure(3);clf;
%    psdplot(S1,'raw');
figure(4);clf;
    %zplot(S0a,S1);   % Compare exact with computed
    zplot(S1);
%figure(5);clf;
%    hplot(S1,[-5,5]);   % Compare exact with computed
