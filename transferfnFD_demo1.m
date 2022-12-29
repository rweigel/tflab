addpath([fileparts(mfilename('fullpath')),'/plot']);
set(0,'DefaultFigureWindowStyle','docked')

%% E(t+1) = B(t)
N = 11;     % Number of input/output points
H = [0,1]'; % Impulse reponse filter E(f) = H(f)*B(f)

% Get in/out signals and exact transfer function
S0 = transferfnFD_demo_signals(-1,struct('H',H,'N',N));

opts = transferfnFD_options(0);         % Use default options
opts.transferfnFD.loglevel = 1;

S1 = transferfnFD(S0.In,S0.Out,opts);   % Compute transfer function

% Add variables to S0 for plotting
%S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
S0.Options.description = 'Actual'; % For plot labels


close all;
% Plot raw input/output data
figure(1)
tsplot(S1,struct('type','raw'));

% Plot actual and predicted output
figure(2)
tsplot(S1,struct('type','error'));

% Plot spectrum of input/output data
figure(3)
psdplot(S1,struct('type','raw'));

% Compare exact Z with computed
figure(4)
zplot({S0,S1});

% Compare exact H with computed
figure(5)
hplot(S1,[-5,5]);   

assert(all(abs(real(S1.Z) - real(S0.Z))) <= 10*eps);
assert(all(abs(imag(S1.Z) - imag(S0.Z))) <= 10*eps);
