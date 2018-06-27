nOrders = 69;
load('ironA.mat') %file containing vector of iron line wavelengths in Angstroms
load('avgWaves.mat')

%making order bins, discretizing into bins
bins = [min(avgWaves, [], 2); max(avgWaves(end, :))];
d = discretize(ironA, bins);

% hasLines = ismember(1:nOrders, d);
coords = [];
for i = 1:nOrders
    orderLines = ironA(d == i);
    for j = 1:length(orderLines)
        [~, col] = min(abs(orderLines(j) - avgWaves(i, :)));
        coords = [coords; i , col];
    end
end
ironAB = ironA(coords(:, 2) > 30 & coords(:,2) < 4066);
blueCoords = coords(coords(:, 2) > 30 & coords(:,2) < 4066, :);

redBins = [min(avgWaves(1, :)); max(avgWaves, [], 2)];
d = discretize(ironA, redBins);

redCoords = [];
for i = 1:nOrders
    orderLines = ironA(d == i);
    for j = 1:length(orderLines)
        [~, col] = min(abs(orderLines(j) - avgWaves(i, :)));
        redCoords = [redCoords; i , col];
    end
end
ironAR = ironA(redCoords(:, 2) > 30 & redCoords(:,2) < 4066);
redCoords = redCoords(redCoords(:, 2) > 30 & redCoords(:,2) < 4066, :);

allCoords = [blueCoords;redCoords];

[coords, iA, ~] = unique(allCoords, 'rows');
ironA = [ironAB; ironAR];
ironA = ironA(iA);

save('ironLineCoords.mat', 'coords', 'ironA');
 