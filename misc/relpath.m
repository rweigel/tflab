function rel = relpath(full)

if length(dbstack) == 0 || isempty(dbstack(1))
    callerpath = pwd;
else
    file = dbstack(1).file;
    callerpath = fileparts(which(file));
end

rel = strrep(full, callerpath, '.');
