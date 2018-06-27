function [serFinal, serErrs] = dateMovingAvg (ser, serE, ts, nDays)
    serFinal = zeros(length(ser), 1);
    serErrs = zeros(length(ser), 1);
    daysInc = zeros(length(ser), 1);
    for i = 1:length(ser)
        thisWindow = ts < ts(i) + days(nDays) & ts >= ts(i);
        [serErrs(i), serFinal(i)] = wmeanSeries(ser(thisWindow), serE(thisWindow));
        daysInc(i) = sum(thisWindow);
    end
end
