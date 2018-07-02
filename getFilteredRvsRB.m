load allFitResults
load dupLines

fSing = f(:, single, :);
fDup = (f(:, pairs(:, 1), :) + f(:, pairs(:, 2), :)) /2;
fComb = [fSing fDup];

offsets = [offsets(single) mean([offsets(pairs(:, 1)); offsets(pairs(:, 2))])];

ironA = [ironA(single); ironA(pairs(:, 1))];
reduced = [reduced(:, single) (reduced(:, pairs(:, 1)) + reduced(:, pairs(:, 2)))/2];

errFit = [errFit(:, single, :) (errFit(:, pairs(:, 1), :) + errFit(:, pairs(:, 2), :))/2];

c = 299792.458 * 1000 ; %speed of light in m/sec

relDepths = squeeze(fComb(:, :, 2)./fComb(:, :, 1));
depthErrs = sqrt(errFit(:, :, 2).^2 + errFit(:, :, 1).^2)./relDepths;
dF = max(relDepths) < 1 & min(relDepths) > 0 & max(depthErrs) < .1; 

gravRed = 636;
rvs = (squeeze(fComb(:, :, 3) + offsets))./ironA' * c - gravRed;
rvErrs = squeeze(errFit(:, :, 3))./ironA'*c;
cF = max(abs(rvs)) < 1000;
sF = std(squeeze(fComb(:, :, 3))./ironA' * c) < 100;
rF = mean(reduced) < 5;

fil =  dF & cF & rF & sF;

rvs = rvs(:, fil);
rvErrs = rvErrs(:, fil);

[~, serPrelim] = wmean(rvs, rvErrs); 
bad = serPrelim - mean(serPrelim) > 10;

rvs = rvs(~bad, :);
rvErrs = rvErrs(~bad, :);

relDepthsG = relDepths(~bad, fil);
depthErrs = depthErrs(~bad, fil);

gIronA = ironA (fil);

save gaussianSeries rvs rvErrs relDepthsG depthErrs gIronA
  
load quadResultsNew

poly3 = @(b, x)(b(1) * x.^3 + b(2) * x.^2 + b(3) * x + b(4));
r = [884.308, - 145.560, - 43.7963, -504.891];
step = .1;
edges = .1:step:.9;

q = [qTot(~bad, single, 3) (qTot(~bad, pairs(:, 1), 3) + qTot(~bad, pairs(:, 2), 3))/2];
polyErr = [errTot(~bad, single, 3) (errTot(~bad, pairs(:, 1), 3) + errTot(~bad, pairs(:, 2), 3))/2];
polyErr = polyErr ./ironA'*c;
polyRVs = q ./ironA' * c - gravRed;

qFil = dF & mean(polyErr) < 50 & max(abs(polyRVs)) < 1000;

relDepthsP = relDepths(:, qFil);
polyRVs = polyRVs(:, qFil);
polyErr = polyErr(:, qFil);
ironA = ironA(qFil);


save rbCorrectedRVs polyRVs polyErr relDepthsP ironA
