function plots(sid,all)

if nargin < 2
    all = 1;
end
addpath([fileparts(mfilename('fullpath')),'/../plot']);

script_dir = fileparts(mfilename('fullpath'));

if ischar(sid)
    cid = sid;
    if strcmp(sid,'KAP')
        cid = 'KAP03'; % TODO: Need to pass cid or find file from dirwalk.
    end
    fname = sprintf('%s/data/%s/%s.mat',script_dir,cid,sid);
    fprintf('Loading %s\n',fname);
    load(fname);
else
    S = sid;
    sid = S{1}.Options.info.stationid;
    cid = S{1}.Options.info.chainid;
end

savedir = sprintf('%s/figures/%s/',script_dir,sid);

opts.type = 'raw';
opts.savefmt = {'pdf'};
opts.filename = [savedir,'timeseries'];
timeseries_plot(S{1},opts);

if all
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
end

% Plot on same axes
opts.filename = [savedir,'SN_compare'];
sn_plot(S,opts);

if all
    opts.savefmt = {'pdf'};
    for i = 1:length(S)
        opts.filename = [savedir,sprintf('transferfnZ_method_%d',i)];
        transferfnZ_plot(S{i},opts);
    end
end

% Plot on same axes
opts.filename = [savedir,'transferfnZ_compare'];
transferfnZ_plot(S,opts);

