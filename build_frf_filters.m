function precomp = build_frf_filters(params)
% BUILD_FRF_FILTERS Precomputes the isotropic spatial-frequency filter bank
% using Will Harrison's polarDistFun.m and createBandPassFilters_customSizes.m.
%
% Converts cpd cutoffs to cycles/image (cpi = cpd * imSize_deg) and optionally
% zeros the DC component of the low-pass band.

imSize      = params.stim.imSize_px;
imSize_deg  = params.stim.imSize_deg;
nyquist_cpd = (imSize / 2) / imSize_deg;

assert(params.sf.highPassFreq_cpd < nyquist_cpd, ...
    'highPassFreq_cpd (%.2f) must be below Nyquist (%.2f cpd).', ...
    params.sf.highPassFreq_cpd, nyquist_cpd);

[~, radDist] = polarDistFun(imSize);

lp_cpi = params.sf.lowPassFreq_cpd  * imSize_deg;
hp_cpi = params.sf.highPassFreq_cpd * imSize_deg;
[sfFilts, steps_cpi] = createBandPassFilters_customSizes( ...
    radDist, params.sf.nFreqs, lp_cpi, hp_cpi);

if params.sf.zeroDC_lp
    lp = sfFilts(:, :, 1);
    lp(radDist == 0) = 0;
    sfFilts(:, :, 1) = lp;
end

sfFiltsShift = zeros(size(sfFilts));
for k = 1:size(sfFilts, 3)
    sfFiltsShift(:, :, k) = fftshift(sfFilts(:, :, k));
end

precomp.radDist      = radDist;
precomp.radDist_cpd  = radDist / imSize_deg;
precomp.sfFilts      = sfFilts;
precomp.sfFiltsShift = sfFiltsShift;
precomp.sfFreqs_cpi  = steps_cpi;
precomp.sfFreqs_cpd  = steps_cpi / imSize_deg;
end
