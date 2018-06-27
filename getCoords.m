load labels
aC = [];
for sC = 1:length(stars)
%     load originalList
%     lineID = 1 + (sC - 1)*nSplit:min([nLines, sC*nSplit]);
%     iA = ironA(lineID);
    load(strcat('lineChunk', num2str(sC), '.mat'))
    ironA = lineChunk;
    save ironA ironA
    findLineCoords
    aC= [aC; coords];
end

save allLineCoords aC