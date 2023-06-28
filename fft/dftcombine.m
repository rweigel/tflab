function DFTc = dftcombine(DFT)
%DFTCOMBINE - Collapse 3-D DFT cell array along 3rd dimension.

% 3rd dim corresponds to segments

DFTc = struct();

for i = 1:size(DFT.In,1)
    tmp = squeeze(DFT.In(i,1,:));
    DFTc.In{i,1} = cat(1,tmp{:});
end

for i = 1:size(DFT.Out,1)
    tmp = squeeze(DFT.Out(i,1,:));
    DFTc.Out{i,1} = cat(1,tmp{:});

    tmp = squeeze(DFT.f(i,1,:));
    DFTc.f{i,1} = cat(1,tmp{:});
end

DFTc.f = DFTc.f{:,1};
DFTc.fe = DFT.fe;