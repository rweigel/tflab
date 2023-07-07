function sdir = scriptdir()
%SCRIPTDIR - Directory of calling script or pwd() if called on command line
%
%   SCRIPTDIR() 

stack = dbstack;

if length(stack) > 1
    sdir = fileparts(which(stack(2).file));
else
    sdir = pwd();
end

