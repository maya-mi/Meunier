load labels
bisectSpanAll = [];
bX = [];
bY = [];
sTot = [];
eTot = [];
for sC = 1:length(stars)
    load(strcat(stars{sC}, '/bisect', stars{sC}))
    bisectSpanAll = [bisectSpanAll bisectSpan];
    bX = cat(2, bX, bisectX);
    bY = cat(2, bY, bisectY);
    sTot = [sTot spanErr];
    eTot = cat(2, eTot, err);
end

save allBisect bisectSpanAll bX bY sTot eTot