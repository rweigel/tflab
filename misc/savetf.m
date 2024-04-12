function savetf(tf,fname)
%SAVETF - Save output of tflab.
%
%   SAVETF(tf,filename) saves the output of tflab to a MATLAB
%   binary file named filename.
%
%   Use tf = load(filename) to recover state.
%
%   See also TFLAB.

if ~exist('fname','var')
    fname = 'tf';
end

base = fileparts(fname);
if ~isempty(base) && ~exist(base, 'dir')
    logmsg('Creating: %s\n', base);
    mkdir(base);
end

fnamefull = fname;
if ~endsWith(fnamefull,'.mat')
    fnamefull = append(fnamefull,'.mat');
end

if isfield(tf,'Segment')
    keeps = {'Z','fe','Metrics'};
    fns = fieldnames(tf.Segment);
    for i = 1:length(fns)
        if ~any(strcmp(fns{i},keeps))
            tf.Segment = rmfield(tf.Segment,fns{i});
        end
    end
end

logmsg('Saving: %s\n',fnamefull);
save(fname,'-v7.3','-struct','tf');
logmsg('Saved:  %s\n',fnamefull);
return

keeps = {'In','Out','Z','dZ','ZVAR','fe','Metrics','Options','Metadata','Segment'};
for i = 1:length(keeps)
    if ~isfield(tf,keeps{i})
        keeps{i} = '';
    end
end

save(fname,'-v7.3','-struct','tf','In','Out',keeps{:});
logmsg('Saved:  %s\n',fnamefull);

