function [B,E,t,outfile,infiles] = EarthScope_read(id)

% List of sites where data are available at baseurl.
prepared = {'VAQ58'};

if all(strcmp(id,prepared) == 0)
    error('Data for %s is not available',id);
end

baseurl = 'http://mag.gmu.edu/git-data/IRIS-EM-download';

outdir  = fullfile(scriptdir(),'data','EarthScope',id);
outfile = sprintf('%s_raw.mat',id);
outfile = fullfile(outdir,outfile);

if strcmp(id,'VAQ58')

    if exist(outfile,'file')
        logmsg('Reading: %s\n',relpath(outfile));
        load(outfile);
        logmsg('Read:    %s\n',relpath(outfile));
        return;
    end

    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    suffix = '2016-06-10_through_2016-06-25.mat';
    channels = {'LFE','LFN','LFZ','LQE','LQN'};
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
        data{c} = load(infile);
    end
    
    B = [data{1}.data,data{2}.data,data{3}.data];
    E = [data{4}.data,data{5}.data];
    t = data{1}.time;

    outfileo = outfile;
    outfile = relpath(outfile);
    logmsg('Saving: %s\n',outfile);
    save(outfileo,'E','B','t','outfile','infiles');
    logmsg('Saved: %s\n',outfile);

end