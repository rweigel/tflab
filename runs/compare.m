% TODO: Add option to show LEMIMT if it was run

scriptpath = fileparts(mfilename('fullpath'));

addpath([fileparts(mfilename('fullpath')),'/../plot']);

fprintf('Reading run files\n');
Middelpos = load([scriptpath,sprintf('/data/Middelpos/Middelpos.mat')]);
KAP103 = load([scriptpath,sprintf('/data/KAP03/KAP103.mat')]);

script_dir = fileparts(mfilename('fullpath'));

savedir = sprintf('%s/figures/KAP103_Middelpos/',script_dir);

popts = struct();

popts.savefmt = {'pdf'};
popts.filename = [savedir,'SN_compare'];
sn_plot({Middelpos.S{1}, Middelpos.S{2}, KAP103.S{3}},popts);

popts.plottype = 1;
popts.filename = [savedir,'transferfnZ_compare'];
transferfnZ_plot({Middelpos.S{1}, Middelpos.S{2}, KAP103.S{3}},popts);
