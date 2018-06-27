load labels
for sC = 1 : length(stars)
    starName = stars{sC};
    load(strcat(starName, '/fitResults', starName, '.mat'))
    load(strcat(starName, '/processed', starName, '.mat'))
    load(strcat(starName, '/propError', starName, '.mat'))
    load(strcat(starName, '/times', starName, '.mat'))

    load(strcat(starName, '/widthsOffsets', starName, '.mat'))

    [nObs, nLines, nPix] = size(normOrders);
    % finalDays = length(edges);
    %Use custom fit for gaussian + vertical shift
    qfit = @(b, x)(b(1) + b(2)*(x - b(3)).^2);
    slopeFit = @(b, x)(b(1) + b(2)*(x - b(4)).^2 + b(3)*(x - b(4)));

    %Pre-designating fit vars
    nNights = length(uniqueNights);
    q = zeros(nNights, nLines, 3);
    e = zeros(nNights, nLines, 3); 
    r = zeros(nNights, nLines);
    m = zeros(nNights, nLines, 2);

    qSlope = zeros(nNights, nLines, 4);
    eQ = zeros(nNights, nLines, 4);
    rQ = zeros(nNights, nLines);
    mQ = zeros(nNights, nLines, 2);


    for i = 1:length(uniqueNights)
        for j = 1:nLines
            try
                thisNight = uniqueNights(i) == obsNights;
                fitX = squeeze(wavelengths(thisNight, j, :)) - ironA(j); 
                [~, lI] = min(abs(mean(fitX) + widths(j)*sqrt(2)/2));
                [~, rI] = min(abs(mean(fitX) - widths(j)*sqrt(2)/2));
                if isempty(lI)
                    lI = 1;
                end
                if isempty(rI)
                    rI = 31;
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

                %figure; plot(fitX, fitY);

                [mn, loc] = min(fitY);
                b0 = [mn, 1, 0];
                [q(i, j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, qfit, b0, 'Weights', (1 ./ errSq));
                e(i, j, :) = sqrt(diag(cov));
                model = qfit(q(i, j, :), fitX);
                chisq = sum((fitY - model) .^2 ./ errSq);
                r(i, j) = chisq ./ (length(fitX) - length(b0));
                qfn = @(x)(squeeze(q(i, j, 1)) + squeeze(q(i, j, 2))*(x - squeeze(q(i, j, 3))).^2);
                [x, fval] = fminbnd(qfn, min(fitX), max(fitX));
                m(i, j, :) = [x, fval];
                

                b0Q = [mn, 1, 0, 0];
                [qSlope(i, j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, slopeFit, b0Q, 'Weights', (1 ./ errSq));
                eQ(i, j, :) = sqrt(diag(cov));
                model = qfit(q(i, j, :), fitX);
                chisqQ = sum((fitY - model) .^2 ./ errSq);
                qfnS = @(x)(squeeze(qSlope(i, j, 1)) + squeeze(qSlope(i, j, 2))*(x - squeeze(qSlope(i, j, 4))).^2 + squeeze(qSlope(i, j, 3))*(x - squeeze(qSlope(i, j, 4))));
                rQ(i, j) = chisqQ ./ (length(fitX) - length(b0Q));
                [x, fval] = fminbnd(qfnS, min(fitX), max(fitX));
                mQ(i, j, :) = [x, fval];

            catch
                q(i, j, :) = nan;
                e(i, j, :) = nan;
                r(i, j) = nan;
                qSlope(i, j, :) = nan;
                eQ(i, j, :) = nan;
                rQ(i, j) = nan;

            end


        end

    end


    save(strcat(starName, '/quadFit', starName, '.mat'), 'q', 'e', 'r', 'qSlope', 'eQ', 'rQ');
    %save('fitResultsHD.mat', 'f', 'errFit', 'chisq', 'reduced', 'ironA')
    clearvars -except starCounter stars
end
