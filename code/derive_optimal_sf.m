function opt = derive_optimal_sf(params)
% DERIVE_OPTIMAL_SF Ideal-observer Fourier analysis to derive the single
% optimal spatial frequency (cpd) for discriminating Karami et al.
% numerosities [6, 10, 17, 29] under a standard human photopic CSF.
%
% For each numerosity, generates balanced displays across all 8 (ItemArea x
% FieldDiam) combinations, computes mean 2D Fourier power |FFT|^2, weights by
% |CSF|^2 (createCSF.m), and finds the peak of |P(N_j) - P(N_i)| for all
% SF > 0.5 cpd across all 6 pairwise numerosity differences. The mean of
% these 6 peak SFs is returned as opt.optimalSF_cpd and saved to
% results/optimal_sf.mat.

if nargin < 1 || isempty(params)
    params = frf_default_params();
end
addpath(params.filterDir);
if ~exist(params.outputDir, 'dir'), mkdir(params.outputDir); end

tStart = tic;
nums   = params.stim.numerosities;
areas  = params.stim.itemAreas_vd2;
diams  = params.stim.fieldDiams_deg;
nNum   = numel(nums);
nPer   = params.csf.nDisplaysPerNum;

% 8 (ItemArea x FieldDiam) combinations
subConds = zeros(numel(areas) * numel(diams), 2);
idx = 0;
for iS = 1:numel(areas)
    for iF = 1:numel(diams)
        idx = idx + 1;
        subConds(idx, :) = [areas(iS), diams(iF)];
    end
end
nSub = size(subConds, 1);

imSize = params.stim.imSize_px;
imDeg  = params.stim.imSize_deg;
[~, radDist] = polarDistFun(imSize);
radCPD = radDist / imDeg;
nyqCPD = (imSize / 2) / imDeg;

% Standard human photopic CSF (createCSF.m)
csf2D = createCSF(max(radCPD, 1e-4), params.csf.fmax, params.csf.vmax, params.csf.beta, params.csf.trunc);
csf2D(radCPD == 0) = 0;
csfPow2D = csf2D .^ 2;

maxBin    = floor(imSize / 2);
fBins_cpd = (1:maxBin) / imDeg;
radRound  = round(radDist);

rng(params.stim.rngSeed);
fprintf('=== Ideal-Observer Optimal SF Derivation (CSF peak=%.1f cpd, SF > %.2f cpd) ===\n', ...
    params.csf.fmax, params.csf.minSF_cpd);

radPow_csf = zeros(nNum, maxBin);
for iN = 1:nNum
    powSum = zeros(imSize, imSize);
    for k = 1:nPer
        si = mod(k - 1, nSub) + 1;
        im = generate_karami_stimulus(nums(iN), subConds(si, 1), subConds(si, 2), params.stim);
        imZero = im - mean(im(:));
        powSum = powSum + abs(fftshift(fft2(imZero))).^2;
    end
    meanPow_csf = (powSum / nPer) .* csfPow2D;
    for b = 1:maxBin
        radPow_csf(iN, b) = mean(meanPow_csf(radRound == b));
    end
    fprintf('  N=%2d (%d displays) done [%.1f s]\n', nums(iN), nPer, toc(tStart));
end

%% Pairwise numerosity differences for SF > minSF_cpd (0.5 cpd)
sfMask     = fBins_cpd > params.csf.minSF_cpd;
fSearch    = fBins_cpd(sfMask);
pairs      = nchoosek(1:nNum, 2);
nPairs     = size(pairs, 1);
pairPeaks  = zeros(nPairs, 1);
pairCurves = zeros(nPairs, maxBin);
pairLabels = cell(nPairs, 1);

fprintf('\nPairwise CSF-weighted Fourier power differences (SF > %.2f cpd):\n', params.csf.minSF_cpd);
for p = 1:nPairs
    i = pairs(p, 1); j = pairs(p, 2);
    diffCurve = abs(radPow_csf(j, :) - radPow_csf(i, :));
    pairCurves(p, :) = diffCurve;
    [~, pkIdx] = max(diffCurve(sfMask));
    pairPeaks(p) = fSearch(pkIdx);
    pairLabels{p} = sprintf('N=%d vs %d', nums(i), nums(j));
    fprintf('  %-12s | peak diff = %.3f cpd\n', pairLabels{p}, pairPeaks(p));
end

optimalSF_cpd = mean(pairPeaks);
fprintf('==> Ideal-observer mean optimal SF = %.3f cpd (SD = %.3f cpd)\n', ...
    optimalSF_cpd, std(pairPeaks));

opt.optimalSF_cpd = optimalSF_cpd;
opt.pairPeaks_cpd = pairPeaks;
opt.pairLabels    = pairLabels;
opt.fBins_cpd     = fBins_cpd;
opt.radPow_csf    = radPow_csf;
opt.pairCurves    = pairCurves;
opt.csfParams     = params.csf;
opt.stimParams    = params.stim;

save(fullfile(params.outputDir, 'optimal_sf.mat'), 'opt');

%% Diagnostic figure
if params.plotSanityChecks
    fig = figure('Name', 'Ideal-Observer Optimal SF Derivation', 'Color', 'w', 'Position', [80 120 1240 460]);
    sfTicks = [0.1, 0.25, 0.5, 1, 2, 4, 8, 16];
    sfTickLabels = {'0.1', '0.25', '0.5', '1', '2', '4', '8', '16'};
    
    % Left: CSF-weighted mean radial power per numerosity
    subplot(1, 2, 1); hold on;
    cols = lines(nNum);
    for iN = 1:nNum
        plot(fBins_cpd, radPow_csf(iN, :), 'Color', cols(iN, :), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('N = %d', nums(iN)));
    end
    xline(params.csf.minSF_cpd, 'k:', '0.5 cpd cutoff', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(optimalSF_cpd, 'r--', sprintf('Optimal: %.2f cpd', optimalSF_cpd), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', sfTicks, 'XTickLabel', sfTickLabels);
    xlim([fBins_cpd(1), nyqCPD]);
    xlabel('Spatial Frequency (cpd)'); ylabel('CSF-Weighted Power |FFT|^2');
    title(sprintf('CSF-Weighted Power Spectra (f_{max}=%.1f cpd)', params.csf.fmax));
    legend('Location', 'southwest'); grid on;
    
    % Right: All 6 pairwise difference curves (> 0.5 cpd) + mean optimal SF
    subplot(1, 2, 2); hold on;
    pCols = parula(nPairs + 1);
    for p = 1:nPairs
        plot(fSearch, pairCurves(p, sfMask), 'Color', pCols(p, :), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%s (%.2f cpd)', pairLabels{p}, pairPeaks(p)));
    end
    xline(optimalSF_cpd, 'r--', sprintf('Mean Optimal = %.2f cpd', optimalSF_cpd), ...
        'LineWidth', 2.0, 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', sfTicks, 'XTickLabel', sfTickLabels);
    xlim([params.csf.minSF_cpd, nyqCPD]);
    xlabel('Spatial Frequency (cpd, > 0.5 cpd)'); ylabel('|\Delta Power|');
    title('Pairwise Numerosity Differences (Ideal-Observer Target SF)');
    legend('Location', 'southwest'); grid on;
    
    outPng = fullfile(params.outputDir, 'fig0_ideal_observer_sf.png');
    if exist('exportgraphics', 'file')
        exportgraphics(fig, outPng, 'Resolution', 150);
    else
        saveas(fig, outPng);
    end
end
end
