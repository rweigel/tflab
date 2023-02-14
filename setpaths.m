function setpaths()
%SETPATHS - Set all paths for transferfnFD
%
%   SETPATHS() Sets all paths

% TODO: Is there a standard way to do this?

setpaths_dir = fileparts(mfilename('fullpath'));

addpath(fullfile(setpaths_dir,'window'));
addpath(fullfile(setpaths_dir,'stats'));
addpath(fullfile(setpaths_dir,'spectra'));
addpath(fullfile(setpaths_dir,'regression'));
addpath(fullfile(setpaths_dir,'plot'));
addpath(fullfile(setpaths_dir,'fft'));
addpath(fullfile(setpaths_dir,'misc'));
addpath(fullfile(setpaths_dir,'lib'));
addpath(fullfile(setpaths_dir,'demos'));
addpath(fullfile(setpaths_dir,'deps/printstruct'));
addpath(fullfile(setpaths_dir,'deps/export_fig'));
