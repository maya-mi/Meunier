starName = 'lineChunk5';

load(strcat(starName, '/fitResults', starName, '.mat'))
load(strcat(starName, '/processed', starName, '.mat'))
load(strcat(starName, '/propError', starName, '.mat'))
load(strcat(starName, '/widthsOffsets', starName, '.mat'))
load(strcat(starName, '/times', starName, '.mat'))


ironAP = ironA;
ironAG = ironA + offsets';

[nObs, nLines, nPix] = size(normOrders);
% finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/(sqrt(2)*b(4))).^2 * .5));
qfit = @(b, x)(b(1) + b(2)*(x - b(3)).^2);

gBeta = zeros(nLines, 4);
pBeta = zeros(nLines, 3);
errG = zeros(nLines, 4);
errPoly = zeros(nLines, 3);
redG = zeros(nLines, 1);
redPoly = zeros(nLines, 1);
m = zeros(nLines, 2);


set(0, 'defaultAxesFontSize', 16)
set(0, 'defaultLineLineWidth', 4)

i = 50;

    for j = 1:nLines
            ironA = ironAG;
            thisNight = uniqueNights(i) == obsNights;
            %Getting full line prof
            allX = squeeze(wavelengths(thisNight, j, :) - ironA(j));
            allX = reshape(allX, 1, numel(allX));
            
            allY = squeeze(normOrders(thisNight, j, :));
            allY = reshape(allY, 1, numel(allY));
            %For Gaussian fit
            fitX = squeeze(wavelengths(thisNight, j, :)) - ironA(j); 
            [~, lI] = min(abs(mean(fitX)  + 1.75 * widths(j) * sqrt(2)));
            [~, rI] = min(abs(mean(fitX)  - 1.75 * widths(j) * sqrt(2)));
            if isempty(lI)
                lI = 1;
            end
            if isempty(rI)
                rI = 1;
            end
            
            fitX = fitX(:, lI:rI);
            fitX = reshape(fitX, 1, numel(fitX));
            [fitX, I] = sort(fitX); %sorting to make bisectors easier
            fitY = squeeze(normOrders(thisNight, j, lI:rI)); 
            fitY = reshape(fitY, 1, numel(fitY));
            fitY = fitY(I);
            err = squeeze(propError(thisNight, j, lI:rI));
            err = reshape(err, 1, numel(err));
            errSq = err(I) .^2;

            [mn, loc] = min(fitY);
            b0 = [1, 1 - mn, fitX(loc), widths(j)];
            [gBeta(j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, mfit, b0, 'Weights', (1 ./ errSq));
            errG(j, :) = sqrt(diag(cov));
            modelG = mfit(gBeta(j, :), fitX);
            chisq = sum((fitY - modelG) .^2 ./ errSq);
            redG(j) = chisq ./ (length(fitX) - length(b0));
            
            %For Polynomial
            ironA = ironAP; 
            fitXP = squeeze(wavelengths(thisNight, j, :)) - ironA(j); 
            [~, lI] = min(abs(mean(fitXP) + widths(j)*sqrt(2)/2));
            [~, rI] = min(abs(mean(fitXP) - widths(j)*sqrt(2)/2));
            if isempty(lI)
                lI = 1;
            end
            if isempty(rI)
                rI = 31;
            end

            fitXP = fitXP(:, lI:rI);
            fitXP = reshape(fitXP, 1, numel(fitXP));
            [fitXP, I] = sort(fitXP); %sorting to make bisectors easier
            fitYP = squeeze(normOrders(thisNight, j, lI:rI)); 
            fitYP = reshape(fitYP, 1, numel(fitYP));
            fitYP = fitYP(I);
            err = squeeze(propError(thisNight, j, lI:rI));
            err = reshape(err, 1, numel(err));
            errSq = err(I) .^2;

            [mn, loc] = min(fitYP);
            b0 = [mn, 1, 0];
            [pBeta(j, :), ~, ~, cov, ~, ~] = nlinfit(fitXP, fitYP, qfit, b0, 'Weights', (1 ./ errSq));
            errPoly(j, :) = sqrt(diag(cov));
            modelP = qfit(pBeta(j, :), fitXP);
            chisq = sum((fitYP - modelP) .^2 ./ errSq);
            redPoly(j) = chisq ./ (length(fitXP) - length(b0));
            qfn = @(x)(squeeze(pBeta(j, 1)) + squeeze(pBeta(j, 2))*(x - squeeze(pBeta(j, 3))).^2);
            [x, fval] = fminbnd(qfn, min(fitXP), max(fitXP));
            m(j, :) = [x, fval];
            
            figure; 
            plot(allX + offsets(j), allY, 'k.', 'MarkerSize', 30)
            hold on
            plot(fitX + offsets(j), modelG, 'r')
            plot(fitXP, modelP, 'b', 'LineWidth', 4.5)
            ylim([0 1])
            xlim([-.2 .2])
            xlabel('Offset from theoretical line center (Angstroms)')
            legend('Observed data point', 'Included in Gaussian fit', 'Included in polynomial fit', 'Location', 'southwest')
            ylabel('Relative Depth')
            title(num2str(ironA(j)))
            saveTitle = sprintf('forNailah/%u_%.0f.jpg', j, ironA(j));
            saveas(gcf, saveTitle)
            close
    end
     
    names = {'Wavelength', 'gaussianFit', 'polyFit', 'gErr', 'pErr', 'gRed', 'pRed', 'pMin'};
    fitParams = table(ironA, gBeta, pBeta, errG, errPoly, redG, redPoly, m, 'VariableNames', names);
    writetable(fitParams, 'forNailah/fitParams.xlsx');
    
    justLambda = table(ironA, 'VariableNames', {'Wavelength'});
    writetable(justLambda, 'forNailah/wavelengths.xlsx')
    
    namesS = {'Wavelength', 'gErr', 'pErr', 'gRed', 'pRed'};
    simp = table(ironA, errG, errPoly, redG, redPoly, 'VariableNames', namesS);
    writetable(simp, 'forNailah/errParams.xlsx')
