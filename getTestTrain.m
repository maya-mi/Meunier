function [trainX, trainY, testX, testY, I] = getTestTrain(X, Y)
    r = rand(length(Y), 1);
    [~, I] = sort(r);

    cutoff = floor(.75*length(Y));

    trainY = Y(I(1:cutoff));
    testY = Y(I(cutoff + 1:end));
    trainX = X(I(1:cutoff), :);
    testX = X(I(cutoff + 1 :end), :);
end