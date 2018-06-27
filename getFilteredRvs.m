load allFitResults
load newSDO

c = 299792.458 * 1000 ; %speed of light in m/sec

dF = squeeze(f(:, :, 2)./ f(:, :, 1));
dF = max(dF) < 1 & min(dF) > 0;

cF = max(abs(squeeze(f(:, :, 3)) + offsets)) < .02;
sF = std(squeeze(f(:, :, 3))./ironA' * c) < 100;
rF = mean(reduced) < 5;

filter =  dF & cF & rF & sF;

gravRed = 636;
rvs = (squeeze(f(:, filter, 3) + offsets(filter)))./ironA(filter)' * c - gravRed;
rvErrs = squeeze(errFit(:, filter, 3))./ironA(filter)'*c;
[nObs, nLines] = size(rvs);
relDepths = squeeze(f(:, filter, 2)./f(:, filter, 1));
depthErrs = sqrt(errFit(:, filter, 2).^2 + errFit(:, filter, 1).^2)./f(:, filter, 1);
widths = squeeze(f(:, filter, 4));

blue = ironA(filter);

[~, serPrelim] = wmean(rvs, rvErrs);
noBadDays = abs(serPrelim - mean(serPrelim)) < 10;

rvsB = rvs(noBadDays, :);
rvErrsB = rvErrs(noBadDays, :);
relDepthsB = relDepths(noBadDays, :);
depthErrsB = depthErrs(noBadDays, :);

% load lineChunk1/timesLineChunk1;
% ts = uniqueNights(noBadDays, :);
% 
% vCon = vCon(noBadDays);

%save blueFits blue rvsB rvErrsB
