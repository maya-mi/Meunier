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

load allLineCoords
coordsB = aC(pairs(:, 1), :);
coordsR = aC(pairs(:, 2), :);

load quadResultsNew
load allBisect
load lineChunk1/timeslineChunk1

c = 299792.458 * 1000 ; %speed of light in m/sec

offB = offsets(pairs(:, 1));
offR = offsets(pairs(:, 2));

% bXB = bX(:, pairs(:, 1), :);
% bYB = bY(:, pairs(:, 1), :);
% bErrB = eTot(:, pairs(:, 1), :);
% 
% bXR = bX(:, pairs(:, 2), :);
% bYR = bY(:, pairs(:, 2), :);
% bErrR = eTot(:, pairs(:, 2), :);

ironAPairs = ironA(pairs(:, 1));

rvsGB = (squeeze(f(:, pairs(:, 1), 3)) + offsets(pairs(:, 1)))./ironAPairs' * c;
rvErrsGB = squeeze(errFit(:, pairs(:, 1), 3))./ironAPairs'*c;

rvsGR = (squeeze(f(:, pairs(:, 2), 3)) + offsets(pairs(:, 2)))./ironAPairs' * c;
rvErrsGR = squeeze(errFit(:, pairs(:, 2), 3))./ironAPairs'*c;

rvsQB = squeeze(qTot(:, pairs(:, 1), 3))./ironAPairs'*c;
rvErrsQB = squeeze(errTot(:, pairs(:, 1), 3))./ironAPairs'*c;

rvsQR = squeeze(qTot(:, pairs(:, 2), 3))./ironAPairs'*c;
rvErrsQR = squeeze(errTot(:, pairs(:, 2), 3))./ironAPairs'*c;

dB = squeeze(f(:, pairs(:, 1), 2)./f(:, pairs(:, 1), 1));
dErrsB = sqrt(errFit(:, pairs(:, 1), 2).^2 + errFit(:, pairs(:, 1), 1).^2).*dB;

dR = squeeze(f(:, pairs(:, 2), 2)./f(:, pairs(:, 2), 1));
dErrsR = sqrt(errFit(:, pairs(:, 2), 2).^2 + errFit(:, pairs(:, 2), 1).^2).*dR;

rR = reduced(:, pairs(:, 1));
rB = reduced(:, pairs(:, 2));

gF = max(abs(rvsGB)) < 1000 & max(abs(rvsGR)) < 1000 & max(rvErrsGB) < 100 & max(rvErrsGR) < 100;
qF = max(abs(rvsQB)) < 1000 & max(abs(rvsQR)) < 1000 & max(rvErrsQB) < 100 & max(rvErrsQR) < 100;
dF = max(dB) < 1 & min(dB) > 0 & max(dR) < 1 & min(dB) > 0;
dEF = max(dErrsB) < .1 & max(dErrsR) < .1;
rF = max(rR) < 5 & max(rB) < 5;
fil = gF & qF & dF & dEF & rF;

iA = ironAPairs(fil);

offB = offB(fil);
offR = offR(fil);

% bXB = bXB(:, fil, :);
% bYB = bYB(:, fil, :);
% bErrB = bErrB(:, fil, :);
% 
% bXR = bXR(:, fil, :);
% bYR = bYR(:, fil, :);
% bErrR = bErrR(:, fil, :);
% 
% save overlapBisect bXB bYB bErrB bXR bYR bErrR offB offR

cB = coordsB(fil, :);
cR = coordsR(fil, :);

rvsQB = rvsQB(:, fil);
rvsQR = rvsQR(:, fil);

rvsGB = rvsGB(:, fil);
rvsGR = rvsGR(:, fil);

rvErrsQB = rvErrsQB(:, fil);
rvErrsQR = rvErrsQR(:, fil);

rvErrsGB = rvErrsGB(:, fil);
rvErrsGR = rvErrsGR(:, fil);

dB = dB(:, fil);
dErrsB = dErrsB(:, fil);

dR = dR(:, fil);
dErrsR = dErrsR(:, fil);

rR = rR(:, fil);
rB = rB(:, fil);

[gSerErrR, gSerR] = wmean(rvsGR, rvErrsGR);
[gSerErrB, gSerB] = wmean(rvsGB, rvErrsGB);

[qSerErrR, qSerR] = wmean(rvsQR, rvErrsQR);
[qSerErrB, qSerB] = wmean(rvsQB, rvErrsQB);

[dSerErrR, dSerR] = wmean(dR, dErrsR);
[dSerErrB, dSerB] = wmean(dB, dErrsB);

save overlap iA cB cR rvsQB rvsQR rvsGB rvsGR rvErrsQB rvErrsQR rvErrsGB rvErrsGR dB dErrsB dR dErrsR rB rR
save overlapSeries gSerErrR gSerR gSerErrB gSerB qSerErrB qSerB qSerErrR qSerR dSerErrR dSerR dSerErrB dSerB
