function [data,outfile,infiles] = EarthScope_read(id)

% List of sites where data are available at baseurl.
prepared = {'VAQ58','ORF03','ORG03'};

if all(strcmp(id,prepared) == 0)
    error('Data for %s is not available',id);
end

baseurl = 'http://mag.gmu.edu/git-data/IRIS-EM-download';

outdir  = fullfile(scriptdir(),'..','data','EarthScope',id);
outfile = sprintf('%s_raw.mat',id);
outfile = fullfile(outdir,outfile);

if exist(outfile,'file')
    logmsg('Reading: %s\n',relpath(outfile));
    load(outfile);
    logmsg('Read:    %s\n',relpath(outfile));
    return;
end

if ~exist(outdir,'dir')
    mkdir(outdir);
end

if strcmp(id,'VAQ58')
    % TODO: Get this from a JSON file posted at
    % http://mag.gmu.edu/git-data/IRIS-EM-download
    % (need to create JSON file with list of sites where data has been
    % downloaded along with their start/stop and channels).
    suffix = '2016-06-10_through_2016-06-25.mat';
    channels = {'LFE','LFN','LFZ','LQE','LQN'};
end

if strcmp(id,'ORF03')
    % TODO: Get this from a JSON file posted at
    % http://mag.gmu.edu/git-data/IRIS-EM-download
    % (need to create JSON file with list of sites where data has been
    % downloaded along with their start/stop and channels).
    suffix = '2007-08-19_through_2007-09-07.mat';
    channels = {'LFE','LFN','LFZ','LQE','LQN'};
end

if strcmp(id,'ORG03')
    % TODO: Get this from a JSON file posted at
    % http://mag.gmu.edu/git-data/IRIS-EM-download
    % (need to create JSON file with list of sites where data has been
    % downloaded along with their start/stop and channels).
    suffix = '2007-08-17_through_2007-09-07.mat';
    channels = {'LFE','LFN','LFZ','LQE','LQN'};
end

[data,infiles] = readFiles(id, channels, suffix);

logmsg('Saving: %s\n',relpath(outfile));
save(outfile,'data','outfile','infiles');
logmsg('Saved: %s\n',relpath(outfile));

function [data,infiles] = readFiles(id, channels, suffix)
    for c = 1:length(channels)
        fname = sprintf('%s-%s-%s',id,channels{c},suffix);
        infile = fullfile(outdir,'measurements',fname);
        infiles{c} = relpath(infile);
        if ~exist(fileparts(infile),'dir')
            mkdir(fileparts(infile));
        end
        url = sprintf('%s/%s/%s',baseurl,id,fname);
        if ~exist(infile,'file')
            logmsg(sprintf('Requesting: %s\n',url));
            websave(infile,url);
            logmsg(sprintf('Wrote: %s\n',relpath(infile)));
        end
        logmsg('Reading: %s\n',relpath(infile));
        data{c} = load(infile);
    end
end
end