load drs
load lineChunk1/timeslineChunk1
load realData

drsRVs = zeros(length(uniqueNights), 1);
drsTimes = datetime(obsJDcut, 'ConvertFrom', 'juliandate');
drsDays = dateshift(drsTimes, 'start', 'day');
err = zeros(length(uniqueNights), 1);
for i = 1:length(uniqueNights)
    try
    thisDay  = uniqueNights(i) == drsDays;
    drsRVs(i) = mean(sunCorrCutDiff(thisDay));
    err(i) = sunRefRVErrCut(thisDay);
    catch
    end
end

fil = ~isnan(drsRVs) & abs(drsRVs - nanmean(drsRVs)) < 8;
t = uniqueNights(fil);
RVdrs = drsRVs(fil);
RVdrs = RVdrs - min(RVdrs);


load newSdoSer

RVconv = RVconv(fil);
RVsppl = RVsppl(fil);

linComb = @(b, x)(b(1) + b(2).*x(:, 1) + b(3).*x(:, 2));
[AB, resid, ~, covFit, ~, ~] = nlinfit([RVconv RVsppl], RVdrs, linComb, [0 1 1]);

RV0 = RV0cen(fil);
[AB0, resid0, ~, covFit0, ~, ~] = nlinfit([RVconv RVsppl], RV0, linComb, [0 1 1]);

linFit = @(b, x)(b(1) * x + b(2));
[rBeta, rRes, ~, rCov, ~, ~] = nlinfit(resid, resid0, linFit, [1 0]);