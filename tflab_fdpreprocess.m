function S = tflab_fdpreprocess(S)

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

[S.DFT.In, S.DFT.f, S.DFT.fe] = dftbands(S.In, opts);
S.DFT.Out = dftbands(S.Out, opts);

all = 1;
if all && isfield(S,'In_')
    fns = fieldnames(S.In_);
    for i = 1:length(fns)
        [dft, f, fe] = dftbands(S.In_.(fns{i}), opts);
        S.DFT.In_.(fns{i}).DFT = dft;
        S.DFT.In_.(fns{i}).f = f;
        S.DFT.In_.(fns{i}).fe = fe;
    end
    fns = fieldnames(S.Out_);
    for i = 1:length(fns)
        [dft, f, fe] = dftbands(S.Out_.(fns{i}), opts);
        S.DFT.Out_.(fns{i}).DFT = dft;
        S.DFT.Out_.(fns{i}).f = f;
        S.DFT.Out_.(fns{i}).fe = fe;
    end
end
