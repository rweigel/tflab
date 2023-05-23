function savetf(tf,fname)
%SAVETF - Save output of tflab.
%
%   SAVETF(tf,filename) saves the output of of tflab to a MATLAB
%   binary file named filename.
%
%   Use tf = load(filename) to recover state.
%
%   See also TFLAB.

base = fileparts(fname);
if ~isempty(base) && ~exist(base, 'dir')
    logmsg(sprintf('Creating: %s\n', base));
    mkdir(base);
end

logmsg(sprintf('Saving: %s.mat\n',fname));

save(fname,'-v7.3','-struct','tf');

logmsg(sprintf('Saved:  %s.mat\n',fname));

