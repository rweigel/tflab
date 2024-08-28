clear;
close all % To reduce memory.
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

print_figs = 0;
Nboot = 100;

for run_num = 0:0
    
    switch run_num
        case 0
            start = '20120712';
            stop = '20120716';
        case 1
            start = '20120712';
            stop = '20120811';
        case 2
            start = '20120712';
            stop = '20121107';
    end
    
    filestr = 'Middelpos';
    dirstr  = sprintf('tfs-%s-%s',start,stop);
    rundir  = fullfile(scriptdir(),'data',filestr,dirstr);
    
    %Middelpos_main(rundir, filestr, start, stop, Nboot);
    Middelpos_plot(rundir, filestr, print_figs);

end
