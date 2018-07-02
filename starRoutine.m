% getTargets
% clear

getCoordBounds

load originalList
nLines = length(ironA);
nSplit = 400;
numSplitFiles = ceil(length(ironA )/nSplit);
stars = cell(1, numSplitFiles);
for sC = 1:numSplitFiles
    lineID = 1 + (sC - 1)*nSplit:min([nLines, sC*nSplit]);
    lineChunk = ironA(lineID);
    save(strcat('lineChunk', num2str(sC), '.mat'), 'lineChunk', 'lineID')
    stars{sC} = strcat('lineChunk', num2str(sC));
end

save labels stars
    
for starCounter = 2:numSplitFiles
    starName = stars{starCounter};
    load targets
    directory = find(strcmp(target, 'Sun'));
    starOffset = 0;
    load(strcat('lineChunk', num2str(starCounter), '.mat'))
    ironA = lineChunk;
    save ironA ironA
    findLineCoords
    mkdir(starName);
%     grabIronLinesSun
%     grabWavesBlazes
    clearvars -except starCounter starName stars rvShift numSplitFiles
    normalizeSun
    clearvars -except starCounter starName stars rvShift numSplitFiles
    firstFit
    if strcmp(starName, 'Sun')
        load Sun/widthsOffsetsSun
        save widths.mat widths
    end
    fit2
    clearvars -except starCounter stars rvShift numSplitFiles
end

grabGaussian
 
quadratic
grabQuads

% bisectors
% grabBisectors
    
postProcess
    
meunierRB   
synthTest
    
    
