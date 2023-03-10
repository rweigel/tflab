function rel = relpath(full)

if length(dbstack) == 0 || isempty(dbstack(1))
    callerpath = pwd;
else
    callerpath = fileparts(which(dbstack(1).file));
end

rel = strrep(full, callerpath, '.');
