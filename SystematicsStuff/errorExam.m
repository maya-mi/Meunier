load allFitResults
c = 299792.458 * 1000 ; %speed of light in m/sec

rvErrs = squeeze(errFit(:, :, 3))./ironA' * c;
depthErrs = squeeze(sqrt(errFit(:, :, 2).^2 + errFit(:, :, 1).^2)./f(:, :, 1));
widthErrs = squeeze(errFit(:, :, 4));

relDepths = squeeze(f(:, :, 2)./f(:, :, 1));

load allCoords

grossF = nanmean(rvErrs) < 50 & nanmean(depthErrs) < .1 & nanmean(widthErrs) < .01;

figure; subplot(2, 3, 1)
plot(coords(grossF, 2), nanmean(depthErrs(:, grossF)), 'o')
title('Pixel Num vs Depth Errs')
ylabel('Normalized Intensity')
xlabel('CCD Pixel Num')
subplot(2, 3, 2)
plot(coords(grossF, 2), nanmean(rvErrs(:, grossF)), 'o')
title('Pixel Num vs RV Errs')
ylabel('RV (m/s)')
xlabel('CCD Pixel Num')
subplot(2, 3, 3)
plot(coords(grossF, 2), nanmean(widthErrs(:, grossF)), 'o')
title('Pixel Num vs Width Errs')
ylabel('Width Error (AU)')
xlabel('CCD Pixel Num')

subplot(2, 3, 4)
plot(coords(grossF, 1), nanmean(depthErrs(:, grossF)), 'o')
title('Order Num vs Depth Errs')
xlim([0 69])
ylabel('Normalized Intensity')
xlabel('Order')
subplot(2, 3, 5)
plot(coords(grossF, 1), nanmean(rvErrs(:, grossF)), 'o')
title('Order Num vs RV Errs')
ylabel('RV (m/s)')
xlabel('Order')
xlim([0 69])
subplot(2, 3, 6)
plot(coords(grossF, 1), nanmean(widthErrs(:, grossF)), 'o')
title('Order Num vs Width Errs')
xlim([0 69])
ylabel('Width Error (AU)')
xlabel('Order')

figure; subplot(1, 3, 1)
plot(ironA(grossF), nanmean(depthErrs(:, grossF)), 'o')
title('Wavelength vs Depth Errs')
ylabel('Normalized Intensity')
xlabel('Angstroms')
xlim([3830 6900])
subplot(1, 3, 2)
plot(ironA(grossF), nanmean(rvErrs(:, grossF)), 'o')
title('Wavelength vs RV Errs')
xlim([3830 6900])
ylabel('RV (m/s)')
xlabel('Angstroms')
subplot(1, 3, 3)
plot(ironA(grossF), nanmean(widthErrs(:, grossF)), 'o')
xlim([3830 6900])
title('Wavelength vs Width Errs')
ylabel('Width Error (AU)')
xlabel('Angstroms')

figure; subplot(1, 2, 1)
histogram(coords(grossF, 2))
xlim([1 4096])
subplot(1, 2, 2)
load avgWaves
plot((diff(avgWaves')'./ avgWaves(:, 1:end-1) * c)')