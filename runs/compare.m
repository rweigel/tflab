% TODO: Add option to show LEMIMT if it was run
set(0,'defaultFigureWindowStyle','docked');

scriptpath = fileparts(mfilename('fullpath'));

addpath([fileparts(mfilename('fullpath')),'/../plot']);

fprintf('Reading run files.\n');
Middelpos = load([scriptpath,sprintf('/data/Middelpos/Middelpos.mat')]);
KAP103 = load([scriptpath,sprintf('/data/KAP03/KAP103.mat')]);

script_dir = fileparts(mfilename('fullpath'));

savedir = sprintf('%s/data/figures/KAP103_Middelpos/',script_dir);

popts = struct();

popts.savefmt = {'pdf'};
popts.filename = [savedir,'SN_compare'];
snplot({Middelpos.S{1}, KAP103.S{1}, KAP103.S{3}},popts);

popts.plottype = 1;
popts.filename = [savedir,'transferfnZ_compare'];
zplot({Middelpos.S{1}, KAP103.S{1}, KAP103.S{3}},popts);
