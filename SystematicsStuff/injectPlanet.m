getFilteredRvs

planet = .6*cos(2*pi/100 * juliandate(ts));

rvsInj = rvs + planet;
[~, meanRv] = wmean(rvsInj, rvErrs);
figure; plot(ts, meanRv)

[pxx, fs] = plomb(meanRv, ts, 2e-6);
figure; plot(fs, pxx)


