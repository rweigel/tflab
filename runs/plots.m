function plots(sid)

addpath([fileparts(mfilename('fullpath')),'/../plot']);

script_dir = fileparts(mfilename('fullpath'));

if ischar(sid)
    fname = sprintf('%s/data/%s/%s.mat',script_dir,sid,sid);
    fprintf('Loading %s\n',fname);
    load(fname);
else
    S = sid;
    sid = S{1}.Options.info.stationid;
end

savedir = sprintf('%s/figures/%s/',script_dir,sid);

opts.type = 'raw';
opts.savefmt = {'pdf'};
opts.filename = [savedir,'timeseries'];
timeseries_plot(S{1},opts);

opts.type = 'error';
opts.savefmt = {'png'};
for i = 1:length(S)
    if ~isempty(S{i})
        opts.filename = [savedir,sprintf('timeseries-error_method_%d',i)];
        timeseries_plot(S{i},opts);
    end
end

opts.savefmt = {'pdf'};
for i = 1:length(S)
    opts.filename = [savedir,sprintf('SN_method_%d',i)];
    sn_plot(S{i},opts);
end

% Plot all
opts.filename = [savedir,'SN_compare'];
sn_plot(S,opts);

opts.savefmt = {'pdf'};
for i = 1:length(S)
    opts.filename = [savedir,sprintf('transferfnZ_method_%d',i)];
    transferfnZ_plot(S{i},opts);
end

% Plot all
opts.filename = [savedir,'transferfnZ_compare'];
transferfnZ_plot(S,opts);

