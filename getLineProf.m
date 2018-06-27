load allFitResults;

relDepths = mean(squeeze(f(:, :, 2)./f(:, :, 1)));
dF = relDepths < .17 & relDepths > .13;
rF = mean(reduced) < 2;
cF = abs(mean(squeeze(f(:, :, 3)) + offsets)) < .02;
F = dF & rF & cF;
[~, In] = sort(reduced(F));
lineIndices = find(F);
j = lineIndices(In(6));
% load originalList
% [~, j] = min(abs(ironA - 6173.3));

load allFitResults
testW = mean(f(:, j, 4));

chunkN = ceil(j/400);
starName = strcat('lineChunk', num2str(chunkN));
newJ = mod(j, 400);

load(strcat(starName, '/processed', starName, '.mat'))
load(strcat(starName, '/propError', starName, '.mat'))
load(strcat(starName, '/times', starName, '.mat'))


fitX = squeeze(wavelengths(:, newJ, :)) - ironA(newJ); 
[lB, lI] = min(abs(mean(fitX)  + 2 * testW * sqrt(2)));
[rB, rI] = min(abs(mean(fitX)  - 2 * testW * sqrt(2)));
    if isempty(lI)
        lI = 1;
    end
    if isempty(rI)
        rI = 31;
    end

fitX = fitX(:, lI:rI);
fitX = reshape(fitX, 1, numel(fitX));

fitY = squeeze(normOrders(:, newJ, lI:rI)); 
fitY = reshape(fitY, 1, numel(fitY));

[fitX, ia, ~] = unique(fitX);
fitY = fitY(ia);

grid = min(fitX):.001:max(fitX);
tempGrid = min(fitX)+.0001:.01:max(fitX);
jagged = interp1(fitX, fitY, grid);
downSample = interp1(grid, jagged, tempGrid);
lineProf = interp1(tempGrid, downSample, grid, 'cubic');

load originalList
save(sprintf('line%.0f', ironA(j)), 'grid', 'lineProf', 'j');


