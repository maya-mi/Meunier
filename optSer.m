function [optSers, indAbs, indNeg, indPos] = optSer(x, xErr, target)
    [~, nParams] = size(x);    
    corrs = corr(x, target);
    [~, Iabs] = sort(abs(corrs), 'descend');
    [~, Ineg] = sort(corrs);
    [~, Ipos] = sort(corrs, 'descend');
    cs = zeros(nParams, 3);
    for i = 1:nParams
        [~, serAbs] = wmean(x(:, Iabs(1:i)), xErr(:, Iabs(1:i)));
        [~, serNeg] = wmean(x(:, Ineg(1:i)), xErr(:, Ineg(1:i)));
        [~, serPos] = wmean(x(:, Ipos(1:i)), xErr(:, Ipos(1:i)));
        cs(i, :) = corr([serAbs serNeg serPos], target);
    end
    %figure; imagesc(cs)
    [~, pref] = max(abs(cs));
    indAbs = Iabs(1:pref(1));
    [~, optAbs] = wmean(x(:, indAbs), xErr(:, indAbs));
    indNeg = Ineg(1:pref(2));
    [~, optNeg] = wmean(x(:, indNeg), xErr(:, indNeg));
    indPos = Ipos(1:pref(3));
    [~, optPos] = wmean(x(:, indPos), xErr(:, indPos));
    
    optSers = [optAbs optNeg optPos];
end
