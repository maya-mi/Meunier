load allFitResults
load dupLines
load newSDO

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

gray = [169 169 169]/250;
figure; plot(mean(relDepthsP), mean(polyRVs), '.', 'Color', gray, 'MarkerSize', 10)

binDepths = edges(2:end) - step/2;
binRVs = zeros(length(edges) - 1, 1);
binStds = zeros(length(edges) - 1, 1);

for i = 1:length(edges) - 1
    thisBin = mean(relDepthsP) > edges(i) & mean(relDepthsP) < edges(i + 1);
    binRVs(i) = mean2(polyRVs(:, thisBin));
    binStds(i) = std(mean(polyRVs(:, thisBin)));
end

hold on
errorbar(binDepths, binRVs, binStds, 'ko-', 'LineWidth', 2, 'CapSize', 0, 'MarkerFaceColor', 'k')
plot(edges, poly3(r, edges), 'r')
legend('Time-avgd values per line', 'Avgd per bin', 'Fit from Reiners 2016', 'Location', 'southeast')
ylim([-1000 400])
xlim([0 1])
xlabel('Normalized Line Depth')
ylabel('Absolute Conv. Blueshift (m/s)')
title('Third Signature in Solar Observations of Fe I lines')
set(gca, 'FontSize', 20)


save rbCorrectedRVs polyRVs polyErr relDepthsP
