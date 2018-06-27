function [err, val] = wmeanSeries (xs, errs)
ws = 1./ errs.^2;
val = nansum(xs.* ws)./nansum(ws);
err = 1./sqrt(nansum(ws));
end
