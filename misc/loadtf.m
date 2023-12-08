function tf = loadtf(fname,variable,preprocess,postprocess)
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
%   tf = LOADTF(filename, '', preprocess, postprocess) executes
%     tf = tflab_preprocess(tf)
%   if preprocess = 1 and
%     tf = tflab_postprocess(tf)
%   if postprocess = 1.
%
%   See also TFLAB.

fnamefull = fname;
if ~endsWith(fnamefull,'.mat')
    fnamefull = append(fnamefull,'.mat');
end

logmsg(sprintf('Reading: %s\n',fnamefull));
if exist('variable','var') && ~isempty(variable)
    tf = load(fname,variable);
    tf = tf.(variable);
    logmsg(sprintf('Read:    ''%s'' from %s\n',variable,fnamefull));
else
    tf = load(fname);
    logmsg(sprintf('Read:    %s\n',fnamefull));
end

if nargin > 2 && preprocess
    tf = tflab_preprocess(tf);
end
if nargin > 3 && postprocess
    tf = tflab_metrics(tf);
end

