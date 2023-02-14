function savetf(S,fname)
%SAVETF - Save output of transferfnFD.
%
%   SAVETF(S,filename) saves the output of of transferfnFD to a MATLAB
%   binary file named filename.
%
%   Use S = load(filename) to recover state.
%
%   See also TRANSFERFNFD.

if ~exist(fileparts(fname),'dir')
    mkdir(fileparts(fname))
end

fprintf('Saving: %s\n',fname);
save(fname,'-v7.3','-struct','S');
fprintf('Saved: %s\n',fname);
