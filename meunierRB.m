load rbCorrectedRVs
load lineChunk1/timeslineChunk1
set(0, 'DefaultLineLineWidth', 2);

c = 299792.458 * 1000 ; %speed of light in m/sec
gravRed = 636;
arbShift = 295;

[~, serPrelim] = wmean(polyRVs, polyErr);
noBadDays = serPrelim - mean(serPrelim) < 20;

rvs = polyRVs(noBadDays, :);
rvErrs = polyErr(noBadDays, :);
t= uniqueNights(noBadDays);
% 
% Centering by min range
[~, RVtest] = wmean(rvs, rvErrs);
RVtest = mean(rvs, 2);
[~, lowestDay] = min(RVtest);
nDays = length(RVtest);
sorted = sort(rvs);
perRow = mean(sorted(100:300, :));

% Meunier method
relDepths = relDepthsP(noBadDays, :);
avgD = mean(relDepths);

edges = .1:.1:.9;
nBins = length(edges) - 1;
nLines = size(rvs, 2);
binDepth = zeros(nBins, 1);
binVal = zeros(nBins, 1);
binErr = zeros(nBins, 1);

for i = 1:nBins
    thisBin = avgD > edges(i) & avgD < edges(i + 1);
    binDepth(i) = mean(avgD(thisBin));
    binVal(i) = mean(perRow(thisBin));
    binErr(i) = std(perRow(thisBin));
end

smoothedSig = interp1(binDepth, binVal, avgD);
thirdSig = abs(perRow - smoothedSig) < 2* mean(binErr);
thirdSig = logical(thirdSig);
cenRVs = rvs - perRow;
% fil = std(cenRVs) < 70 & thirdSig;
% cenRVs = cenRVs(:, fil);
% avgD = avgD(fil);
% rvErrs = rvErrs(:, fil);

 

s0 = avgD > .2 & avgD < .9;
s1 =  avgD > .5 & avgD < .9; 
s2 = avgD > .2 & avgD < .5;

graphTitles = "S1 = shallow";

[err0, RV0cen] = wmean(cenRVs(:, s0), rvErrs(:, s0));
[err1, RV1cen] = wmean(cenRVs(:, s1), rvErrs(:, s1));
[err2, RV2cen] = wmean(cenRVs(:, s2), rvErrs(:, s2));


[~, RV0] = wmean(rvs(:, s0), rvErrs(:, s0));
[~, RV1] = wmean(rvs(:, s1), rvErrs(:, s1));
[~, RV2] = wmean(rvs(:, s2), rvErrs(:, s2));

RV0NW = mean(rvs(:, s0), 2);
RV1NW = mean(rvs(:, s1), 2);
RV2NW = mean(rvs(:, s2), 2);

RV0NWcen = mean(cenRVs(:, s0), 2);
RV1NWcen = mean(cenRVs(:, s1), 2);
RV2NWcen = mean(cenRVs(:, s2), 2);
% RV0cen = RV0 - min(RV0);
% RV1cen = RV1 - min(RV1);


RVis = [RV1cen RV2cen];
RVerrs = [err1 err2];
nSubs = size(RVerrs, 2);

alphaGuess = .5:.01:1.5;
alphaCorrs = zeros(length(alphaGuess), 2);
meanRVSup = zeros(length(alphaGuess), 2);
stdRat = zeros(length(alphaGuess), 2);
minStd = zeros(length(alphaGuess), 2);
stdSmoothed = zeros(length(alphaGuess), 2);
smoothedRat = zeros(length(alphaGuess), 2);

for j = 1:nSubs
    RVsub = RVis(:, j);
    err = RVerrs(:, j);

    for i = 1:length(alphaGuess)
        alpha = alphaGuess(i);
        vConvOpt = (RVsub - RV0cen)./(alpha - 1);
        vSupOpt = (alpha * RV0cen - RVsub)/(alpha - 1);
        vSupErr = (err0 *alpha - err)/(alpha - 1);
        %Taking lin comb, not quadrature, bc errs not independent
        alphaCorrs(i, j) = corr(vConvOpt, vSupOpt);
        meanRVSup(i, j) = mean(vSupOpt);
        stdRat(i, j) = std(vConvOpt)/std(vSupOpt);
        minStd(i, j) = std(vSupOpt);
        [supSmooth, ~] = dateMovingAvg(vSupOpt, vSupErr, t, 30);
        stdSmoothed(i, j) = std(supSmooth);
        smoothedRat(i, j) = std(vConvOpt)/std(supSmooth);
    end
end

%Method 1
[~, j] = min(abs(meanRVSup));
alphaMean = alphaGuess(j);
%Method 2
[~, m] = min(abs(alphaCorrs));
alphaCorr = alphaGuess(m);
%Method 3
noOff = @(b, x) b*x;
alphaSlope = [nlinfit(RV0cen, RV1cen, noOff, 1), nlinfit(RV0cen, RV2cen, noOff, 1)];

%Method 4
[~, k] = max(stdRat);
alphaStd = alphaGuess(k);

%Method 5
[~, n] = min(stdSmoothed);
alphaSmooth = alphaGuess(n);
% 

%Method 6?
[~, p] = max(smoothedRat);
alphaSmoothRat = alphaGuess(p);

%Method 7?
alphaAbs = [mean(RV1)/mean(RV0), mean(RV2)/mean(RV0)];


[alphaMean; alphaCorr; alphaSlope; alphaStd; alphaSmooth; alphaSmoothRat; alphaAbs]

set(0,'defaultAxesFontSize',14)
gray = [169 169 169]/250;
orange = [255 147 0]/255;
plot(t, vConvOpt, 'o', 'Color', orange)
hold on
plot(t, vSupOpt, 'o', 'Color', gray)
title('Reconstructed Timeseries: Alpha = .83')
ylabel('RV m/s')
legend('RVconv', 'RVsppl')


save realData err0 err1 err2 RV0 RV1 RV2 RV0cen RV1cen RV2cen t cenRVs rvErrs s0 s1 s2

