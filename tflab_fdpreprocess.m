function S = tflab_fdpreprocess(S)

if isfield(S,'Options')
    opts = S.Options;
else
    opts = tflab_options(0);
end

if isfield(S,'In_')
    In = S.In_.Final;
    Out = S.Out_.Final;
else
    In = S.In;
    Out = S.Out;
end

[S.DFT.In, S.DFT.f, S.fe] = dftsegments(In, opts);
S.DFT.Out = dftsegments(Out, opts);
S.DFT.Weights = dftweights(S.DFT.f, S.DFT.In, S.DFT.Out,opts);
