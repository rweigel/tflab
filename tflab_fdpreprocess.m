function S = tflab_fdpreprocess(S)

opts = S.Options;

[S.DFT.In, S.DFT.f, S.DFT.fe] = dftbands(S.In, opts);
S.DFT.Out = dftbands(S.Out, opts);

all = 1;
if all && isfield(S,'In_')
    fns = fieldnames(S.In_);
    for i = 1:length(fns)
        [S.DFT.In_.(fns{i}),S.DFT.f,S.DFT.fe] = dftbands(S.In_.(fns{i}), opts);
    end
end
if all && isfield(S,'Out_')    
    fns = fieldnames(S.Out_);
    for i = 1:length(fns)
        S.DFT.Out_.(fns{i}).DFT = dftbands(S.Out_.(fns{i}), opts);
    end
end
