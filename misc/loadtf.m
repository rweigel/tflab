function tf = loadtf(fname,variable)
%LOADTF   Load output of tflab.
%
%   tf = LOADTF(filename) load all output of a tflab run.
%
%   var = LOADTF(filename, varstring) load a top-level structure variable.
%
%   Example:
%       In  = LOADTF(filename, 'In');
%       Out = LOADTF(filename, 'Out');
%
%   See also TFLAB.

fnamefull = fname;
if ~endsWith(fnamefull,'.mat')
    fnamefull = append(fnamefull,'.mat');
end

logmsg(sprintf('Reading: %s\n',fnamefull));
if exist('variable','var')
    tf = load(fname,variable);
    tf = tf.(variable);
    logmsg(sprintf('Read:    ''%s'' from %s\n',variable,fnamefull));
else
    tf = load(fname);
    logmsg(sprintf('Read:    %s\n',fnamefull));
end

