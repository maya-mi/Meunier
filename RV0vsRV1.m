load realData

set(0,'defaultAxesFontSize',14)

[absR, absErr, absRSlope,absErrSlope] = makePlot(RV0, RV1, 'Absolute RVs', 'RV1 saturated (m/s)', 'abs1.jpg');

[R, err, RSlope,errSlope] = makePlot(RV0cen, RV1cen, 'RVi(t) - minVal(i)', 'RV1 saturated (m/s)', 'cen1.jpg');

makePlot(RV0, RV2, 'Absolute RVs', 'RV2 shallow (m/s)', 'abs2.jpg');

makePlot(RV0cen, RV2cen, 'RVi(t) - minVal(i)', 'RV2 shallow (m/s)', 'cen2.jpg');


function [R, err, RSlope,errSlope] = makePlot(RV0, RV1, t, ylab, saveT)
figure;
gray = [169 169 169]/250;
orange = [255 147 0]/255;
plot(RV0, RV1, 'o', 'Color', gray)
lin = @(b, x)(b(1)*x + b(2));
[beta,R,~,CovB,~,~] = nlinfit(RV0, RV1, lin, [1 0]);
err = sqrt(diag(CovB));
hold on
plot(RV0, RV0*beta(1) + beta(2), 'k')
title(t)
xlabel('RV0 (m/s)')
ylabel(ylab)
axis equal
fitLabel = sprintf('RV1 = RV0 * %.2f + %.2f', beta(1), beta(2));
%text(max(RV0) - 6, min(RV1) +1, fitLabel, 'FontSize', 20)
noOff = @(b, x) b*x;
[betaSlope,RSlope,~,CovBSlope,~,~] = nlinfit(RV0, RV1, noOff, 1);
errSlope = sqrt(diag(CovBSlope));
fitLabelSlope = sprintf('RV1 = RV0 * %.2f', betaSlope);
plot(RV0, RV0*betaSlope, 'Color', orange)
legend('Values per time', fitLabel, fitLabelSlope, 'Location', 'southeast')
saveas(gcf, saveT)
end