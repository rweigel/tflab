clear;
scriptpath = fileparts(mfilename('fullpath'));

Middelpos = load([scriptpath,sprintf('/data/Middelpos/Middelpos.mat')]);
KAP103 = load([scriptpath,sprintf('/data/KAP03/KAP103.mat')]);

popts = struct();
    popts.savefmt = struct();
    popts.savefmt.pdf = 0;
    popts.savedir = 'figures/KAP103_Middelpos';

timeseries_plot({KAP103.S{4},KAP103.S{1}},popts);

opts.plottype = 1;
%opts.period_range = [4,1e4];
transferfnZ_plot({KAP103.S{4},KAP103.S{1},Middelpos.S{1}},popts);
