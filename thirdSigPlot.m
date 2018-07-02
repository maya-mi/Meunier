load rbCorrectedRVs

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