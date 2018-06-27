load allFitResults
[n, bin] = histc(ironA, unique(ironA));
single = n == 1;
multiple = find(n > 1);
index    = find(ismember(bin, multiple));
pairs = [];
for i = 1:length(ironA)
    pairI = find(ironA == ironA(i));
    if length(pairI)  == 2
        pairs = [pairs; pairI'];
    end  
end
pairs = unique(sort(pairs, 2, 'descend'), 'rows');

save dupLines pairs multiple single