clear;
close all % To reduce memory.
addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

rerun       = 0;
run_nums    = [5];
print_figs  = 1;
%Nboot       = 100;
Nboot       = 0;

for run_num = run_nums

    switch run_num
        case 0 % 4 days
            start = '20120712';
            stop = '20120716';
        case 1 % 1 month
            start = '20120712';
            stop = '20120811';
        case 2 % All data (note significant data gaps; don't use)
            start = '20120712';
            stop = '20121107';
        case 3
            start = '20170101'; % Problem with data after this date
            stop = '20170104';
        case 4 % 4 days
            start = '20120714';
            stop = '20120718';
        case 5 % All days before file with much missing data
            start = '20120712';
            stop = '20120904';
    end
    
    filestr = 'Middelpos';
    dirstr  = sprintf('tfs-%s-%s',start,stop);
    rundir  = fullfile(scriptdir(),'data',filestr,dirstr);
    
    if rerun == 1 || ~exist(rundir, 'dir')
        Middelpos_main(rundir, filestr, start, stop, Nboot);
    end
    Middelpos_plot(rundir, filestr, print_figs);

end
