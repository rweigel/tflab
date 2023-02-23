%% Transfer function is time shift
% E(t+1) = B(t)

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

N = 11;     % Number of input/output points
H = [0,1]'; % Impulse reponse filter E(f) = H(f)*B(f)

% Get in/out signals and exact transfer function
S0 = demo_signals('fromH',struct('H',H,'N',N));

opts = tflab_options(0);         % Use default options
opts.tflab.loglevel = 1;

S1 = tflab(S0.In,S0.Out,opts);   % Compute transfer function

% Add variables to S0 for plotting
%S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0
S0.Options.description = 'Actual'; % For plot labels


close all;
% Plot raw input/output data
figure()
    tsplot(S1,struct('type','raw'));

% Plot actual and predicted output
figure()
    tsplot(S1,struct('type','error'));

% Plot spectrum of input/output data
figure()
psdplot(S1,struct('type','raw'));

% Compare exact Z with computed
figure()
zplot({S0,S1});

% Compare exact H with computed
figure()
hplot(S1,[-5,5]);   

assert(all(abs(real(S1.Z) - real(S0.Z))) <= 10*eps);
assert(all(abs(imag(S1.Z) - imag(S0.Z))) <= 10*eps);
