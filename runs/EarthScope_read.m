function [B,E,t,matfile] = EarthScope_read(id)

baseurl = 'http://mag.gmu.edu/git-data/IRIS-EM-download';

% Dir of data file
outdir = fullfile(scriptdir(),'data',id);

if strcmp(id,'VAQ58')

    matfile = sprintf('%s_raw.mat',id);
    matfile = fullfile(outdir,matfile);
    if exist(matfile,'file')
        logmsg('Reading: %s\n',relpath(matfile));
        load(matfile);
        logmsg('Read:    %s\n',relpath(matfile));
        return;
    end

    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    suffix = '2016-06-10_through_2016-06-25.mat';
    channels = {'LFE','LFN','LFZ','LQE','LQN'};
    for c = 1:length(channels)
        fname = sprintf('%s-%s-%s',id,channels{c},suffix);
        outfile = fullfile(outdir,'measurements',fname);
        if ~exist(fileparts(outfile),'dir')
            mkdir(fileparts(outfile));
        end
        url = sprintf('%s/%s/%s',baseurl,id,fname);
        if ~exist(outfile,'file')
            websave(outfile,url);
        else
            data{c} = load(outfile);
        end
    end
    
    B = [data{1}.data,data{2}.data,data{3}.data];
    E = [data{4}.data,data{5}.data];
    t = data{1}.time;

    logmsg('Saving: %s\n',matfile);
    save(matfile,'E','B','t','matfile');
    logmsg('Saved: %s\n',matfile);

end