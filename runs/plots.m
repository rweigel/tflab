function plots(sid,all)

if nargin < 2
    all = 1;
end
addpath([fileparts(mfilename('fullpath')),'/../plot']);

script_dir = fileparts(mfilename('fullpath'));

if ischar(sid)
    cid = sid;
    if strncmp(sid,'KAP',3)
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

savedir = sprintf('%s/data/figures/%s/',script_dir,sid);

opts.savefmt = {'pdf'};

opts.type = 'raw';
opts.filename = [savedir,'timeseries'];
figure()
tsplot(S{1},opts);

% Plot on same axes
figure()
opts.filename = [savedir,'SN_compare'];
[ax1,ax2] = snplot(S,opts);

% Plot on same axes
opts.filename = [savedir,'transferfnZ_compare'];
figure()
[ax1, ax2] = zplot(S,opts);


if all
    opts.type = 'error';
    opts.savefmt = {'png'};
    for i = 1:length(S)
        if ~isempty(S{i})
            figure()
            opts.filename = [savedir,sprintf('timeseries-error_method_%d',i)];
            tsplot(S{i},opts);
        end
    end

    opts.savefmt = {'pdf'};
    for i = 1:length(S)
        opts.filename = [savedir,sprintf('SN_method_%d',i)];
        figure()
        snplot(S{i},opts);
    end
    
    opts.savefmt = {'pdf'};
    for i = 1:length(S)
        opts.filename = [savedir,sprintf('transferfnZ_method_%d',i)];
        figure()
        zplot(S{i},opts);
    end
end

