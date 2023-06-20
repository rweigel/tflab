function DFTc = dftcombine(DFT)
%DFTCOMBINE

DFTc = struct();
for i = 1:size(DFT.In,1)

    tmp = squeeze(DFT.In(i,1,:));
    DFTc.In{i,1} = cat(1,tmp{:});

    tmp = squeeze(DFT.Out(i,1,:));
    DFTc.Out{i,1} = cat(1,tmp{:});

    tmp = squeeze(DFT.f(i,1,:));
    DFTc.f{i,1} = cat(1,tmp{:});

    tmp = squeeze(DFT.InWeights(i,1,:));
    DFTc.InWeights{i,1} = cat(1,tmp{:});

    tmp = squeeze(DFT.OutWeights(i,1,:));
    DFTc.OutWeights{i,1} = cat(1,tmp{:});
end

DFTc.fe = DFT.fe;