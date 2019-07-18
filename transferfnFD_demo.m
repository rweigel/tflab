clear;

addpath([fileparts(mfilename('fullpath')),'/plot']);
addpath([fileparts(mfilename('fullpath')),'/..']);

close all;
set(0,'defaultFigureWindowStyle','docked');

% Create signals and return H used to compute them along with Z.
S0 = transferfnFD_demo_signals(0,1);
S0.description = 'Actual';

opts = transferfnFD_options(0);
S1 = transferfnFD(S0.In,S0.Out,opts);

S0.Time = S1.Time;                 % Use default computed time from S1 for S0
S0.Options.info = S1.Options.info; % Use same variable labels from S1 for S0

timeseries_plot(S1,'raw'); % Plot raw input/output data
transferfnZ_plot(S0,S1);   % Compare exact with computed
transferfnH_plot(S0,S1);   % Compare exact with computed
