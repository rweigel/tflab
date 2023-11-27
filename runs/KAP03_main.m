% Add location of tflab_setpaths() to path.
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
% Set all tflab paths.
tflab_setpaths(); 

% Get input/output data
[B,E,t,infile,outfile,timedelta] = KAP03_clean(stationid);

% Set variable information (used for plots)
meta = struct();
    meta.instr     = {'$B_x$','$B_y$'};
    meta.inunit    = 'nT';
    meta.outstr    = {'$E_x$','$E_y$'};
    meta.outunit   = 'mV/km';
    meta.timeunit  = 's';
    meta.timedelta = timedelta;
    meta.timestart = datestr(t(1),'yyyy-mm-ddTHH:MM:SS.FFF');
    meta.frequnit  = 'Hz';
    meta.freqsf    = 1/timedelta; % Multiply frequencies by this to get frequnit
    meta.chainid   = chainid;
    meta.stationid = stationid;
