function DFTc = dftcombine(DFT)
%DFTCOMBINE - Collapse 3-D DFT cell array along 3rd dimension.

% 3rd dim corresponds to segments

DFTc = struct();

for i = 1:size(DFT.In,1)
    s = size(DFT.In{i});
    DFTc.In{i,1} = transpose(reshape(DFT.In{i}, s(2), s(1)*s(3)));
end

for i = 1:size(DFT.Out,1)
    s = size(DFT.Out{i});
    DFTc.Out{i,1} = transpose(reshape(DFT.Out{i}, s(2), s(1)*s(3)));
    DFTc.f{i,1} = repmat(DFT.f{i,1}, s(3), 1);
end

DFTc.fe = DFT.fe;