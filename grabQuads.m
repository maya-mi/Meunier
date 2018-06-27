load originalList

qTot = [];
errTot = [];
rTot = [];
load labels

for sC = 1:length(stars)
    starName = stars{sC};
    load(strcat(starName, '/quadFit', starName, '.mat'))
    qTot = cat(2, qTot, q);
    errTot = cat(2, errTot, e);
    rTot = cat(2, rTot, r);
end
    

save quadResultsNew qTot errTot rTot
