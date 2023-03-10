function savetf(TF,fname)
%SAVETF - Save output of tflab.
%
%   SAVETF(TF,filename) saves the output of of tflab to a MATLAB
%   binary file named filename.
%
%   Use TF = load(filename) to recover state.
%
%   See also TFLAB.

if ~exist(fileparts(fname),'dir')
    mkdir(fileparts(fname))
end

logmsg(sprintf('Saving: %s\n',fname));

save(fname,'-v7.3','-struct','TF');

logmsg(sprintf('Saved: %s\n',fname));

