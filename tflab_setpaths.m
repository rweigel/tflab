function tflab_setpaths()
%TFLAB_SETPATHS - Set all paths for tflab
%
%   TFLAB_SETPATHS() Sets all paths for TFLAB.
%
%   See also TFLAB.

% TODO: Determine if there is standard way to set paths for a package.

if 0
    % TODO: Move this to tflab.m and have it check if functions from
    % either toolbox is needed before proceeding.
    if isempty(which('regress'))
      % Trigger missing Statistics Toolbox error message with an install link
      help('regress');
    end
    if isempty(which('rectwin'))
      % Trigger missing Signal Processing Toolbox message with an install link
      help('regress');
    end
end

% TODO?: Check if paths already set and return.

% Directory of tflab_setpaths.
pkgdir = fileparts(mfilename('fullpath'));

addpath(fullfile(pkgdir,'data'));
addpath(fullfile(pkgdir,'demos'));
addpath(fullfile(pkgdir,'deps'));
addpath(fullfile(pkgdir,'deps/printstruct'));
addpath(fullfile(pkgdir,'deps/export_fig'));
addpath(fullfile(pkgdir,'filter'));
addpath(fullfile(pkgdir,'lib'));
addpath(fullfile(pkgdir,'fft'));
addpath(fullfile(pkgdir,'misc'));
addpath(fullfile(pkgdir,'plot'));
addpath(fullfile(pkgdir,'plot/adjust'));
addpath(fullfile(pkgdir,'plot/annotate'));
addpath(fullfile(pkgdir,'regression'));
addpath(fullfile(pkgdir,'spectra'));
addpath(fullfile(pkgdir,'stats'));
addpath(fullfile(pkgdir,'window'));
