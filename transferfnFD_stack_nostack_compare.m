addpath([fileparts(mfilename('fullpath')),'/plot']);
set(0,'DefaultFigureWindowStyle','docked')

%% Test of E(t+1) = B(t)
N = 1000;     % Number of input/output points
H = [0,1]'; % Impulse reponse filter E(f) = H(f)*B(f)

% Get in/out signals and exact transfer function
S0 = transferfnFD_demo_signals(-1,struct('H',H,'N',N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = transferfnFD_options(0); % Default options
opts.transferfnFD.loglevel = 1;
opts.td.window.width = N/10; 
opts.td.window.shift = N/10; 

S1 = transferfnFD(S0.In,S0.Out,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('---\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = transferfnFD_options(0); % Default options
opts.transferfnFD.loglevel = 1;
opts.fd.stack.average.function = '';

S2 = transferfnFD(S0.In,S0.Out,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
close all;
% Plot raw input/output data
tsplot(S1,struct('type','raw'));

% Plot actual and predicted output
tsplot(S1,struct('type','error'));

% Plot spectrum of input/output data
psdplot(S1,struct('type','raw'));

% Compare exact Z with computed
zplot({S0,S1});

% Compare exact H with computed
hplot(S1,[-5,5]);   

assert(all(abs(real(S1.Z) - real(S0.Z))) <= 10*eps);
assert(all(abs(imag(S1.Z) - imag(S0.Z))) <= 10*eps);
end

