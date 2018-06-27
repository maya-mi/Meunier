load allFitResults
load newSDO

c = 299792.458 * 1000 ; %speed of light in m/sec

dF = squeeze(f(:, :, 2)./ f(:, :, 1));
dF = max(dF) < 1 & min(dF) > 0;

cF = max(abs(squeeze(f(:, :, 3)) + offsets)) < .02;
sF = std(squeeze(f(:, :, 3))./ironA' * c) < 100;
rF = mean(reduced) < 5;

filter =  dF & cF & rF & sF;

gravRed = 636;
rvs = (squeeze(f(:, filter, 3) + offsets(filter)))./ironA(filter)' * c - gravRed;
rvErrs = squeeze(errFit(:, filter, 3))./ironA(filter)'*c;
[nObs, nLines] = size(rvs);
relDepths = squeeze(f(:, filter, 2)./f(:, filter, 1));
depthErrs = sqrt(errFit(:, filter, 2).^2 + errFit(:, filter, 1).^2)./f(:, filter, 1);
widths = squeeze(f(:, filter, 4));

ironA = ironA(filter);

%Constructing RV corrs; 
rvCorr = corrcoef(rvs);
rvCorr(rvCorr == 1) = 0; %Set diag elements to 0. 
rvCorr(abs(rvCorr) < .5) = 0; 
RV_G = graph(rvCorr);
bins = conncomp(RV_G);
%figure; plot(RV_G);

%test = subgraph(RV_G, find(bins == mode(bins)));
%figure; plot(test);

tbl = tabulate(bins); 


nonTriv = unique(tbl(tbl(:, 2) ~= 1, 1));
rvGrouped = zeros(nObs, length(nonTriv));
for i = 1:length(nonTriv)
    thisBin = bins == nonTriv(i);
    [~, rvGrouped(:, i)] = wmean(rvs(:, thisBin), rvErrs(:, thisBin)); 
end

cenRVGrouped = rvGrouped - mean(rvGrouped); 

%Doing the same with depths
depthCorr = corrcoef(relDepths);
depthCorr(depthCorr == 1) = 0; %Set diag elements to 0. 
depthCorr(depthCorr < .1) = 0; 
D_G = graph(depthCorr);
bins_d = conncomp(D_G);
%figure; plot(RV_G);

test = subgraph(D_G, find(bins_d == mode(bins_d)));
%figure; plot(test);

tbl_d = tabulate(bins_d); 

nonTriv_d = unique(tbl_d(tbl_d(:, 2) ~= 1, 1));
depthGrouped = zeros(nObs, length(nonTriv_d));
for i = 1:length(nonTriv_d)
    thisBin = bins_d == nonTriv_d(i);
    [~, depthGrouped(:, i)] = wmean(relDepths(:, thisBin), depthErrs(:, thisBin)); 
end

    
load lineChunk1/timeslineChunk1
figure; histogram(corr(depthGrouped, vCon), -1:.1:1)
%[~, test] = wmean(relDepths, depthErrs);
    

X = [depthGrouped, cenRVGrouped];
Y = vCon;
r = rand(length(Y), 1);
[~, I] = sort(r);

cutoff = floor(.75*length(Y));

trainY = Y(I(1:cutoff));
testY = Y(I(cutoff + 1:end));

trainX = X(I(1:cutoff), :);
testX = X(I(cutoff + 1 :end), :);

mdl = TreeBagger(1000, trainX, trainY, 'Method', 'regression', 'OOBPrediction', 'on', 'OOBPredictorImportance', 'on');
testTree = predict(mdl, testX);
figure; scatter(testY, 2* testTree)
figure; bar(mdl.OOBPermutedVarDeltaError)