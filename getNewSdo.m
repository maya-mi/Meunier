load area
load lineChunk1/timeslineChunk1
[nObs, nFeat] = size(area); 
tSDO = datetime(area(:, 1), 'ConvertFrom', 'juliandate');
t = dateshift(tSDO, 'start', 'day');
allFeat = zeros(length(uniqueNights), nFeat - 1);
labels = {'|B|', 'Activity Fraction', 'Activity Fraction (Faculae)', 'Activity Fraction (Spots)', 'Vcon', 'Vphot', 'Vquiet'};
for i = 1:length(uniqueNights)
    allFeat(i, :) = mean(area(t == uniqueNights(i), 2:end));
end

vCon = allFeat(:, 5); 
save newSDO vCon allFeat labels