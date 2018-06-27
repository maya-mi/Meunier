% load lineChunk1/timesLineChunk1
% nDays = length(uniqueNights);
% 
% load newSDO
% 
% SDOdays = dateshift(datetime(jTime, 'ConvertFrom', 'juliandate'), 'start', 'day');
% RVsppl = zeros(nDays, 1);
% RVconvNoArea = zeros(nDays, 1);
% 
% load areaSDO
% 
% SDOdaysArea = dateshift(datetime(jvC, 'ConvertFrom', 'juliandate'), 'start', 'day');
% RVconv = zeros(nDays, 1);
% for i = 1:nDays
%     thisDay = SDOdays == uniqueNights(i);
%     thisDayArea = SDOdaysArea == uniqueNights(i);
%     RVsppl(i) = mean(vSppl(thisDay));
%     RVconv(i) = mean(vC(thisDayArea));
%     RVconvNoArea(i) = mean(vCNA(thisDay));  
% 
% end
% 
% save newSdoSer RVsppl RVconv RVconvNoArea

load newSdoSer
alpha1 = .8;
alpha2 = 1.2;

load realData;

rRV0 = mean(cenRVs(:, s0), 2);
rRV1 = mean(cenRVs(:, s1), 2);
rRV2 = mean(cenRVs(:, s2), 2);


nDays = length(t);
norm0 = abs(random('Normal', 0, median(err0) , nDays, 1));
norm1 = abs(random('Normal', 0, median(err1) , nDays, 1));
norm2 = abs(random('Normal', 0, median(err2) , nDays, 1));


R = corrcoef([err0 err1 err2]);
errMat = zeros(3, 3, nDays);
errMat(1, 1, :) = err0;
errMat(2, 2, :) = err1;
errMat(3, 3, :) = err2;
cov = zeros(3, 3, nDays);
for i = 1:nDays
    cov(:, :, i) = errMat(:, :, i)*R*errMat(:,:, i);
end
ser = mvnrnd([0 0 0], cov);

linComb = @(b, x)(b(1) + b(2).*x(:, 1) + b(3).*x(:, 2));
[AB, resid, ~, covFit, ~, ~] = nlinfit([RVconv RVsppl], RV0cen, linComb, [0 1 1]);

RVconv = RVconv*AB(2);
RVsppl = RVsppl*AB(3);


RV1 = RVsppl + alpha1*RVconv + ser(:, 2);
RV2 = RVsppl + alpha2*RVconv + ser(:, 3); 
RV0 = RVsppl + RVconv + ser(:, 1);
SerErr0 = err0;
SerErr1 = err1;
SerErr2 = err2;

gray = [169 169 169]/250;
orange = [255 147 0]/255;
figure; plot(t, RVconv, 'o', 'Color', orange)
hold on
plot(t, RVsppl, 'o', 'Color', gray)
plot(t, ser(:, 1), 'ko')
legend('RVconv', 'RVsppl', 'b0')
ylabel('RV (m/s)')
saveas(gcf, 'sdo.jpg')
figure; 
plot(t, RV0, 'ko')
ylabel('RV (m/s)')
title('RV0: RVsppl + RVconv + b0(t)')
saveas(gcf, 'rv0Synth.jpg')
figure; 
plot(t, RV0cen, 'ko')
ylabel('RV (m/s)')
title('RV0: real')
saveas(gcf, 'rv0Real.jpg')
% plot(t, RV0cen, 'o')

RVis = [RV1 RV2];
RVerrs = [SerErr1 SerErr2];
nSubs = size(RVerrs, 2);

alphaGuess = .5:.05:1.5;
alphaCorrs = zeros(length(alphaGuess), 2);
meanRVSup = zeros(length(alphaGuess), 2);
stdRat = zeros(length(alphaGuess), 2);
minStd = zeros(length(alphaGuess), 2);
stdSmoothed = zeros(length(alphaGuess), 2);
smoothedRat = zeros(length(alphaGuess), 2);
% 
% for j = 1:nSubs
%     RVsub = RVis(:, j);
%     err = RVerrs(:, j);
% 
%     for i = 1:length(alphaGuess)
%         alpha = alphaGuess(i);
%         vConvOpt = (RVsub - RV0)./(alpha - 1);
%         vSupOpt = (alpha * RV0 - RVsub)/(alpha - 1);
%         vSupErr = (SerErr0 *alpha - err)/(alpha - 1);
%         %Taking lin comb, not quadrature, bc errs not independent
%         alphaCorrs(i, j) = corr(vConvOpt, vSupOpt);
%         meanRVSup(i, j) = mean(vSupOpt);
%         stdRat(i, j) = std(vConvOpt)/std(vSupOpt);
%         minStd(i, j) = std(vSupOpt);
%         [supSmooth, ~] = dateMovingAvg(vSupOpt, vSupErr, t, 30);
%         stdSmoothed(i, j) = std(supSmooth);
%         smoothedRat(i, j) = std(supSmooth)/std(vConvOpt);
%     end
% end
% 
% %Method 1
% [~, j] = min(abs(meanRVSup));
% alphaMean = alphaGuess(j);
% %Method 2
% [~, m] = min(abs(alphaCorrs));
% alphaCorr = alphaGuess(m);
% %Method 3
% alphaSlope = [RV0\RV1, RV0\RV2];
% 
% %Method 4
% [~, k] = max(stdRat);
% alphaStd = alphaGuess(k);
% 
% %Method 5
% [~, n] = min(stdSmoothed);
% alphaSmooth = alphaGuess(n);
% % 
% 
% %Method 6?
% [~, p] = min(smoothedRat);
% alphaSmoothRat = alphaGuess(p);
% 
% 
% [alphaMean; alphaCorr; alphaSlope; alphaStd; alphaSmooth; alphaSmoothRat]



function [err] = makeSynthScatter(fitErr)
    nDays = length(fitErr);
    err = zeros(nDays, 1);
    for i = 1:nDays
        err(i) = random('Normal', 0, fitErr(i));
    end
end

