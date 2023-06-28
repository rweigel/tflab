function tf = loadtf(fname,variable)
%SAVETF - Load output of tflab.
%
%   tf = LOADTF(filename) load the output.
%
%   See also TFLAB.

fnamefull = fname;
if ~endsWith(fnamefull,'.mat')
    fnamefull = append(fnamefull,'.mat');
end

logmsg(sprintf('Reading: %s\n',fnamefull));
if exist('variable','var')
    tf = load(fname,variable);
else
    tf = load(fname);
end
logmsg(sprintf('Read:    %s\n',fnamefull));

