load originalList

% linSamp = linspace(min(ironA), max(ironA), length(ironA));
% binEdges = linspace(min(ironA), max(ironA), length(ironA)/100);
% %figure; histogram(linSamp, binEdges)
% ironA = linSamp;
% findLineCoords;
% linCoords = coords;
% save linTest linSamp linCoords

load allCoords
load avgWaves
nOrders = 69;

set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',20)

starts = avgWaves(:, 1);
ends = avgWaves(:, end);
overlap = ends(1:end - 1) - starts(2:end);
overlap = [overlap; 0];
figure; 
subplot(1, 2,2)
yyaxis left
histogram(coords(:, 1))
ylabel('Number of Lines')
yyaxis right
plot(ends - starts - overlap)
ylabel('Wavelength covered in order (no overlap)')
ylim([0, 84])
xlabel('Order')
xlim([1 69])
title('Order distribution, Fe I')

subplot(1, 2, 1)
histogram(ironA, linspace(min(ironA), max(ironA), length(ironA)/200))
xlim([min(ironA), max(ironA)])
ylabel('Number of Lines')
xlabel('Wavelength (Angstroms)')
title('Wavelength Distribution, Fe I')

pixelOver = zeros(nOrders, 1);
for i = 1:nOrders- 1
    [~, pixelOver(i)] = min(abs(avgWaves(i, :) - starts(i + 1)));
end
pixelOver(end) = 4096;

figure; 
subplot(1, 2, 1)
plot(coords(:, 2), coords(:, 1), '.')
set (gca,'YDir','reverse')
hold on
plot(pixelOver, 1:1:69)
xlim([31 4065])
xlabel('Pixel Num')
ylim([1 69])
ylabel('Order Num')
lgd = legend('Line Coord','Order Cut-off (Pixel Num)', 'Location', 'southwest');
lgd.FontSize = 16;
title('Order/Pixel Coordinates, Fe I')
subplot(1,2 ,2 )
histogram(coords(:, 2), linspace(31, 4065, 15))
xlim([31 4065])
xlabel('Pixel Num')
ylabel('Number of lines')
title('Distribution across pixels, Fe I')

figure; plot(ends - starts, 'k:')
hold on
plot(overlap, 'k--')
plot(ends - starts - overlap)
xlabel('Order')
xlim([1 69])
legend('Total wavelength range', 'Overlap from following order', 'Wavelength covered, no overlap', 'Location', 'northwest')
title('Wavelength Covered Per Order')
ylabel('Angstroms')
