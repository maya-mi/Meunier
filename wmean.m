function [serErrs, ser] = wmean(xs, errs)
ws = 1./ nanmean(errs).^2;
ser = nansum(xs.* ws, 2)./nansum(ws);
serErrs = sqrt(nansum(errs.^2.*ws.^2, 2))./nansum(ws);
end

