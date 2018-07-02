load overlap
load overlapSeries
load lineChunk1/timeslineChunk1
set(0,'defaultAxesFontSize',14)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor'},{'k','k'})

%Making figs
% figure; histogram(cB(:, 2), 'FaceColor', 'b')
% hold on
% histogram(cR(:, 2), 'FaceColor', 'r')
% title('Pixel Numbers of Overlapped Lines')
% xlabel('Pixel')
% ylabel('Num Lines')
% xlim([0, 4069])
% saveas(gcf, 'hist.jpg')

param = {'d', 'gc', 'poly'};

bs = {dB, rvsGB, rvsQB};
bEs = {dErrsB, rvErrsGB, rvErrsQB};

rs = {dR, rvsGR, rvsQR};
rEs = {dErrsR, rvErrsGR, rvErrsQR};

axLab = {'Rel. Depth (Normalized Intensity)', 'RV (m/s)', 'RV(m/s)'};
tHist = {'Relative Depth', 'Gaussian Line Centers', '2nd-degree poly line centers'};
tPix = cellfun(@(s) ['Pixel Dep. of ' s], tHist, 'UniformOutput', false);
tPlot = cellfun(@(s) ['Error-weighted ' s], tHist, 'UniformOutput', false);
tOrd = cellfun(@(s) ['Order Dep. of ' s], tHist, 'UniformOutput', false);

bad = abs(gSerB - mean(gSerB)) > 10;
t = juliandate(uniqueNights);

rLim = [min(cR(:, 2)) 4096];
bLim = [0, max(cB(:, 2))];

for i = 3:length(param)
%     makePlot(bs{i}, rs{i}, bEs{i}, rEs{i}, t, bad, tPlot{i})
%     
 pixDep(bs{i}, rs{i}, bEs{i}, rEs{i}, bLim, rLim, cB, cR, axLab{i}, tPix{i})
%     
%     orderDep(bs{i}, rs{i}, bEs{i}, rEs{i}, cB, cR, axLab{i}, tOrd{i})
    
   % makeHist(bs{i}, rs{i}, bEs{i}, rEs{i}, axLab{i}, tHist{i}, bad)
%     set(groot,{'DefaultAxesXColor','DefaultAxesYColor'},{'k','k'})
end



function [] = makePlot(b, r, bErr, rErr, times, bad, tPlot)
    [bSerErr, bSer] = wmean(b, bErr);
    [rSerErr, rSer] = wmean(r, rErr);
    figure; errorbar(times, bSer, bSerErr, 'b.')
    hold on
    errorbar(times, rSer, rSerErr, 'r.')
    blueLab = sprintf('Blue fit: %.3g std', std(bSer(~bad)));
    redLab = sprintf('Red fit: %.3g std', std(rSer(~bad)));
    legend(blueLab, redLab)
    title(tPlot)
end

function [R] = makeHist(b, r, bErr, rErr, axisLab, tHist, bad)
    b = b(~bad, :);
    r = r(~bad, :);
    bErr = bErr(~bad, :);
    rErr = rErr(~bad, :);
    
    [bSerErr, bSer] = wmean(b, bErr);
    [rSerErr, rSer] = wmean(r, rErr);
    figure; 
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor'},{'b','r'})
    subplot(1, 3, 1)
    plot(b, r, '.', 'MarkerEdgeColor', [169,169,169]/250)
    hold on
    mB = mean(b);
    mR = mean(r);
    
    ws = 1./ (std(b).^2 + std(r).^2);
    line = @(b, x) (b(1) + b(2) *x);
    [beta, R, ~, cov, ~, ~] = nlinfit(mB, mR, line, [0 1], 'Weights', ws);
    err = sqrt(diag(cov));
    plot(mB, beta(1) + beta(2) * mB, 'k--')
    plot(mB, mR, 'ko', 'MarkerFaceColor', 'k')
    
    title(sprintf('Time-avgd val per line: %.2f (+/- %.3f) + %.2f (+/- %.2f)', beta(2), err(2), beta(1), err(1)))
    
    xlabel(axisLab)
    ylabel(axisLab)

    subplot(1, 3, 2)
    errorbar(bSer, rSer, rSerErr, rSerErr, bSerErr, bSerErr, 'k.')
    xlabel(axisLab)
    ylabel(axisLab)
    title('Error-weighted mean timeseries')
   

    subplot(1, 3, 3)
    nLines = size(bErr, 2);
    meanErrsB = sqrt(sum(bErr.^2, 2))/nLines; %Adding in quadrature
    meanErrsR = sqrt(sum(rErr.^2, 2))/nLines;
    errorbar(mean(b, 2), mean(r, 2), meanErrsR, meanErrsR, meanErrsB, meanErrsB, 'k.')
    xlabel(axisLab)
    ylabel(axisLab)
    title('Unweighted mean timeseries')
    
    x = suptitle(tHist);
    set(x, 'FontSize', 20, 'FontWeight', 'bold')
end


function [] = pixDep(b, r, bE, rE, bLim, rLim, cB, cR, ylab, tPix)
    figure; 
    subplot(2, 2, 1)
    plot(cB(:, 2), mean(bE), 'bo')
    xlim(bLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Time averaged mean error per line, blue')

    subplot(2, 2, 2)
    plot(cR(:, 2), mean(rE), 'ro')
    xlim(rLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Time averaged mean error per line, red')

    subplot(2, 2, 3)
    plot(cB(:, 2), std(b), 'bo')
    xlim(bLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Timeseries std per line, blue')

    subplot(2, 2, 4)
    plot(cR(:, 2), std(r), 'ro')
    xlim(rLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Timeseries std per line, red')
    
    suptitle(tPix)
    
    figure; 
    subplot(1,2,1)
    plot(cB(:, 2), std(b) ./ mean(bE) , 'bo')
    ylabel('Ratio (unitless)')
    xlabel('Pixel Number')
    xlim(bLim)
    title('Timeseries std / time-averaged fit error ratio')
    subplot(1, 2, 2)
    plot(cR(:, 2), std(r) ./ mean(rE) , 'ro')
    ylabel('Ratio (unitless)')
    xlabel('Pixel Number')
    xlim(rLim)
    title('Timeseries std / time-averaged fit error ratio')
    suptitle(tPix)
    
    figure; 
    subplot(1, 2, 1)
    plot(cB(:, 2), std(b), 'bo')
    xlim(bLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Timeseries std per line, blue')

    subplot(1, 2, 2)
    plot(cR(:, 2), std(r), 'ro')
    xlim(rLim)
    ylabel(ylab)
    xlabel('Pixel Number')
    title('Timeseries std per line, red')
end

function [] = orderDep(b, r, bE, rE, cB, cR, ylab, tOrd)
figure; 
subplot(1, 2, 1)
plot(cB(:, 1), mean(bE), 'bo')
hold on
plot(cR(:, 1), mean(rE), 'ro')
xlim([0 69])
ylabel(ylab)
title('Mean error per line')

subplot(1, 2, 2)
plot(cB(:, 1), sqrt(mean((r - b).^2)), 'ko')
xlim([0 69])
ylabel(ylab)
title('RMS(redVal(t) - blueVal(t)) per line')

x = suptitle(tOrd);
set(x, 'FontSize', 20, 'FontWeight', 'bold')
end