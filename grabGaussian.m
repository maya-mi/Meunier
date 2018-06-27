load originalList
nLines = length(ironA);
nSplit = 400;
numSplitFiles = ceil(length(ironA )/nSplit);
fCombined = [];
wCombined = [];
oCombined = [];
iACombined = [];
r = [];
eF = [];
totFilt = [];
iATest = [];
for sC = 1:numSplitFiles
    load(strcat('lineChunk', num2str(sC), '.mat'))
    ironA = lineChunk;
    save ironA ironA
    findLineCoords
    iA = ironA;
    load(strcat(stars{sC}, '/fitResults', stars{sC}))
    load(strcat(stars{sC}, '/widthsOffsets', stars{sC}))
    fCombined = cat(2, fCombined, f);
    wCombined = [wCombined widths];
    oCombined = [oCombined offsets];
    iACombined = [iACombined; iA];
    iATest = [iATest; ironA];
    r = cat(2, r, reduced);
    eF = cat(2, eF, errFit);
end

f = fCombined;
widths = wCombined;
offsets = oCombined;
ironA = iACombined;
reduced = r;
errFit = eF;

save allFitResults f widths offsets ironA reduced errFit totFilt