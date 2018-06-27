load allBisect
load originalList
load lineChunk1/timeslineChunk1
c = 299792.458 * 1000 ; %speed of light in m/sec

load allFitResults
relDepths = mean(squeeze(f(:, :, 2)./f(:, :, 1)));
rF = mean(reduced) < 2;
cF = abs(mean(squeeze(f(:, :, 3)) + offsets)) < .02;
edges = .1:.2:.9;
titles = {'.1 < Rel. Depth < .3', '.3 < Rel. Depth < .5', '.5 < Rel. Depth < .7', '.7 < Rel. Depth < .9'};
for depthBin = 4:length(edges) - 1
    dF = edges(depthBin) < relDepths & relDepths < edges(depthBin + 1);
    F = dF & rF & cF;
    % [~, In] = sort(reduced(F));
    % lineIndices = find(F);
    % j = lineIndices(In(6));
    [~, j] = min(abs(ironA - 6173.3));

    X = bX(:, F, :);
    Y = bY(:, F, :);
    err = eTot(:, F, :);
    [~, ser] = wmean(X, err);

    x = squeeze(ser);
    y = squeeze(nanmean(Y));

    % spans = bisectSpanAll(:, j);
    % spanErr = sTot(:, j);

    %load line6173

    % figure; plot(grid, lineProf,'k', 'LineWidth', 3)
    % hold on
    % plot(nanmean(x), nanmean(y), 'k--')
    % title('Line Profile: 6173.3 A')
    % xlabel('Wavelength (angstroms)')
    % ylabel('Normalized intensity')

    %figure; errorbar(x(1, :), y(1, :), err(1, :))

    [coeff,score,latent,tsquared,explained,mu] = pca(x);

    plotX = x./mean(ironA(F)) * c;
    figure;
    for pC = 1:4
        subplot(2, 2, pC)
        plotNum(plotX, y, coeff, explained, pC)
    end
    suptitle(titles{depthBin})
end

%     %Training to vCon
%     load lineChunk1/sdoComp
%     vCon = vs(:, 1);
% 
%     load DRS_RVs_Cut_FWHMCorr
%     tSDO = datetime(obsJDcut, 'ConvertFrom', 'juliandate');
%     t = dateshift(tSDO, 'start', 'day');
%     ccfB = zeros(length(uniqueNights), 1);
% 
%     for i = 1:length(uniqueNights)
%         ccfB(i) = mean(bis_spancut(t == uniqueNights(i)));
%     end
% 
% 
%     r = rand(length(ccfB), 1);
%     [~, I] = sort(r);
% 
%     cutoff = floor(.75*length(ccfB));
%     trainX = x(I(1:cutoff), :);
%     trainY = ccfB(I(1:cutoff));
%     testX = x(I(cutoff + 1:end), :);
%     testY = ccfB(I(cutoff + 1:end));
%     trainXAlt = score(I(1:cutoff), 1:4);
%     testXAlt = score(I(cutoff+1:end), 1:4);
%     %Linear regression model
%     mdl = fitlm(trainX, trainY);
%     pred = predict(mdl, testX);
% 
%     tMdl = TreeBagger(1000, trainX, trainY, 'Method', 'regression', 'OOBPrediction', 'on', 'OOBPredictorImportance', 'on');
%     tPred = predict(tMdl, testX);
%     %figure; bar(tMdl.OOBPermutedVarDeltaError)
% 
%     tMdlAlt = TreeBagger(1000, trainXAlt, trainY, 'Method', 'regression', 'OOBPrediction', 'on', 'OOBPredictorImportance', 'on');
%     tPredAlt = predict(tMdlAlt, testXAlt);
%     %figure; bar(tMdlAlt.OOBPermutedVarDeltaError)
    
function [] = plotNum(x, y, coeff, explained, num)
   set(0, 'defaultAxesFontSize', 20)
   col = repmat([0 0 1], length(coeff), 1);
   col(coeff(:, num) < 0, :) = repmat([1 0 0], sum(coeff(:, num) < 0), 1);
   scatter(nanmean(x), nanmean(y), abs(coeff(:, num)*500), col, 'filled')
   xlim([-50 150])
   xlabel('RV (m/s)')
   ylabel('Normalized Intensity')
   titleSpec = 'PC %d: %.2f Explained';
   title(sprintf(titleSpec, num, explained(num)));
end
