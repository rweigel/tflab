
break





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

tsplot(S1,'raw'); % Plot raw input/output data
zplot(S0,S1);   % Compare exact with computed
hplot(S0,S1,[-10,10]);   % Compare exact with computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
