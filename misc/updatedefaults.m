function opts = updatedefaults(opts,iopts)
%UPDATEDEFAULTS Recursively update select values in a structure
%
%  s = UPDATEDEFAULTS(s, snew) If s is the default structure, then field
%  elements in snew are applied to s.
%
%  Example: s    = struct('field1',11,'field2',22);
%           snew = struct('field1',-11);
%
%           updatedefaults(s, snew) is struct('field1',-11,'field2',22);
%
%

fns = fieldnames(iopts);
for i = 1:length(fns)
    if isfield(opts,fns{i})
        if isstruct(iopts.(fns{i}))
            opts.(fns{i}) = updatedefaults(opts.(fns{i}),iopts.(fns{i}));
        else
            opts.(fns{i}) = iopts.(fns{i});
        end
    else
        warning('Option %s is not a valid option',fns{i});
    end
end
