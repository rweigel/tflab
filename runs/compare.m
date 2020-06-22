%clear;
close all;
set(0,'defaultFigureWindowStyle','docked');
scriptpath = fileparts(mfilename('fullpath'));

addpath([fileparts(mfilename('fullpath')),'/../']);
addpath([fileparts(mfilename('fullpath')),'/../window']); 
addpath([fileparts(mfilename('fullpath')),'/../fft']); 
addpath([fileparts(mfilename('fullpath')),'/../lib']); 
addpath([fileparts(mfilename('fullpath')),'/../plot']);

Middelpos = load([scriptpath,sprintf('/data/Middelpos/Middelpos.mat')]);
KAP103 = load([scriptpath,sprintf('/data/KAP03/KAP103.mat')]);

popts = struct();
    popts.savefmt = struct();
    popts.savefmt.pdf = 0;
    popts.savedir = 'figures/KAP103_Middelpos';

for i = 1:length(KAP103.S)
    KAP103.S{i}.Options.info.stationid = 'KAP103';
end
for i = 1:length(Middelpos.S)
    Middelpos.S{i}.Options.info.stationid = 'Middelpos';
end

timeseries_plot({KAP103.S{4}, KAP103.S{1}},popts);

sn_plot({KAP103.S{4}, KAP103.S{1}, KAP103.S{2}},popts);
    
opts.plottype = 1;
%opts.period_range = [4,1e4];
transferfnZ_plot({KAP103.S{4},KAP103.S{1}, KAP103.S{2}},popts);
