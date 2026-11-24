%% PLOT_IDEAL_OBSERVER_WALKTHROUGH
% Generates step-by-step figures illustrating the 5 stages of the
% ideal-observer spatial-frequency model and its integration into MEG RSA:
%   1) Test displays (Karami et al. dot arrays)
%   2) Filter test displays with human photopic CSF model (createCSF.m)
%   3) Convert to frequency domain (2D & radial |FFT|^2 * CSF^2)
%   4) Look for indicative differences in the spectrum (SF > 0.5 cpd -> 2.56 cpd)
%   5) Use those values in the RSA analysis (2.56 cpd RDM & MEG tuning)
%
% Outputs (saved to FRF_model/results/):
%   - fig_ideal_observer_pipeline.png       (5-step horizontal overview)
%   - fig_ideal_observer_steps1to3.png      (N=6 vs N=29: display -> CSF -> FFT)
%   - fig_ideal_observer_steps4to5.png      (2D/1D spectral diff -> RSA RDM & tuning)

clear; close all; clc;

baseDir = fileparts(mfilename('fullpath'));
resDir  = fullfile(baseDir, 'results');
suppDir = fullfile(resDir, 'panels_supplement');
params  = frf_default_params();
addpath(params.filterDir);
if ~exist(resDir, 'dir'), mkdir(resDir); end
if ~exist(suppDir, 'dir'), mkdir(suppDir); end

% Load precomputed ideal-observer & model RDM outputs
optFile = fullfile(resDir, 'optimal_sf.mat');
assert(exist(optFile, 'file') == 2, 'Run derive_optimal_sf.m (or runAll.m) first.');
load(optFile, 'opt');

rdmFile = fullfile(resDir, 'frf_rdms.mat');
assert(exist(rdmFile, 'file') == 2, 'Run run_frf_rdm_simulation.m (or runAll.m) first.');
load(rdmFile, 'rdm', 'sfFreqs_cpd');

sfSumFile = fullfile(resDir, 'sf_band_summary.csv');
hasMegSum = exist(sfSumFile, 'file') == 2;
if hasMegSum
    sfTbl = readtable(sfSumFile);
end

%% Generate two representative test displays (N = 6 vs N = 29, matched AIA & TFA)
rng(params.stim.rngSeed);
exAIA = 0.12; % deg^2
exTFA = 13.5; % deg
im6   = generate_karami_stimulus(6,  exAIA, exTFA, params.stim);
im29  = generate_karami_stimulus(29, exAIA, exTFA, params.stim);

imSize = params.stim.imSize_px;
imDeg  = params.stim.imSize_deg;
[~, radDist] = polarDistFun(imSize);
radCPD = radDist / imDeg;

% 2D photopic CSF (DC = 0)
csf2D = createCSF(max(radCPD, 1e-4), params.csf.fmax, params.csf.vmax, params.csf.beta, params.csf.trunc);
csf2D(radCPD == 0) = 0;

% Step 2: Filter displays with CSF in Fourier domain -> back to spatial domain
F6  = fftshift(fft2(im6  - mean(im6(:))));
F29 = fftshift(fft2(im29 - mean(im29(:))));

im6_csf  = real(ifft2(ifftshift(F6  .* csf2D)));
im29_csf = real(ifft2(ifftshift(F29 .* csf2D)));
csfLim   = max([abs(im6_csf(:)); abs(im29_csf(:))]) * 0.85;

% Step 3: 2D CSF-weighted power spectra
P6_csf  = (abs(F6)  .* csf2D) .^ 2;
P29_csf = (abs(F29) .* csf2D) .^ 2;

% Frequency axis for zoomed 2D plots (+/- 6 cpd)
fAxis = ((-imSize/2):(imSize/2 - 1)) / imDeg;
zLim  = 6.0;
zMask = abs(fAxis) <= zLim;
fZoom = fAxis(zMask);

% Circle coordinates for optimal SF = 2.562 cpd
theta = linspace(0, 2*pi, 360);
optSF = opt.optimalSF_cpd;
cx = optSF * cos(theta);
cy = optSF * sin(theta);

% Center band index matching optimalSF_cpd
[~, optBandIdx] = min(abs(sfFreqs_cpd - optSF));

%% ========================================================================
%  FIGURE 1: 5-Step Pipeline Overview (1 x 5)
%  ========================================================================
fig1 = figure('Visible', 'off', 'Name', 'Ideal Observer 5-Step Pipeline', 'Color', 'w', 'Position', [30 180 1750 380]);

% Step 1: Test Display
subplot(1, 5, 1);
imshow(im29, [0 1]);
title({'1. Test Displays', 'e.g. N=29 (0.12 deg^2, 13.5\circ)'}, 'FontSize', 10);

% Step 2: Filter with CSF
subplot(1, 5, 2);
imagesc(im29_csf, [-csfLim, csfLim]); axis image off; colormap(gca, gray(256));
title({'2. CSF-Filtered Display', sprintf('createCSF (f_{max}=%.1f cpd)', params.csf.fmax)}, 'FontSize', 10);

% Step 3: Frequency Domain (Radial CSF-weighted spectra)
subplot(1, 5, 3); hold on;
cols4 = lines(4);
nums  = params.stim.numerosities;
for iN = 1:numel(nums)
    plot(opt.fBins_cpd, opt.radPow_csf(iN, :), 'Color', cols4(iN, :), 'LineWidth', 1.7, ...
        'DisplayName', sprintf('N=%d', nums(iN)));
end
xline(params.csf.minSF_cpd, 'k:', 'LineWidth', 1.1, 'HandleVisibility', 'off');
xline(optSF, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.25 0.5 1 2 4 8 16], ...
    'XTickLabel', {'0.25','0.5','1','2','4','8','16'});
xlim([opt.fBins_cpd(1), 16]);
xlabel('Spatial Frequency (cpd)'); ylabel('CSF Power |FFT|^2');
title({'3. Frequency Domain', 'Radial CSF-Weighted Power'}, 'FontSize', 10);
legend('Location', 'southwest', 'FontSize', 8); grid on; axis square;

% Step 4: Indicative Spectral Differences
subplot(1, 5, 4); hold on;
sfMask  = opt.fBins_cpd > params.csf.minSF_cpd;
fSearch = opt.fBins_cpd(sfMask);
nPairs  = size(opt.pairCurves, 1);
pCols   = parula(nPairs + 1);
for p = 1:nPairs
    plot(fSearch, opt.pairCurves(p, sfMask), 'Color', pCols(p, :), 'LineWidth', 1.5);
end
xline(optSF, 'r--', sprintf('Mean Peak\n%.2f cpd', optSF), 'LineWidth', 1.8, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.5 1 2 4 8 16], ...
    'XTickLabel', {'0.5','1','2','4','8','16'});
xlim([params.csf.minSF_cpd, 16]);
xlabel('Spatial Frequency (cpd)'); ylabel('|\Delta Power|');
title({'4. Spectral Differences', '6 Pairwise Peaks (>0.5 cpd)'}, 'FontSize', 10);
grid on; axis square;

% Step 5: RSA Analysis (2.56 cpd RDM)
subplot(1, 5, 5);
imagesc(rdm.SF_Bands{optBandIdx}); axis square; colorbar;
title({'5. Use in RSA Analysis', sprintf('%.2f cpd Band RDM (32\\times32)', sfFreqs_cpd(optBandIdx))}, 'FontSize', 10);
xlabel('Condition'); ylabel('Condition');

save_fig_png(fig1, fullfile(suppDir, 'supp_ideal_observer_pipeline.png'));
close(fig1);

%% ========================================================================
%  FIGURE 2: Detailed Steps 1-3 (N=6 vs N=29 Spatial & Fourier Walkthrough)
%  ========================================================================
fig2 = figure('Visible', 'off', 'Name', 'Ideal Observer Steps 1-3 Detail', 'Color', 'w', 'Position', [50 90 1520 720]);

% Top row: N = 6
subplot(2, 4, 1);
imshow(im6, [0 1]);
title('1a. Test Display (N = 6)', 'FontSize', 11);

subplot(2, 4, 2);
f1D = logspace(log10(0.05), log10(16), 300);
csf1D = createCSF(f1D, params.csf.fmax, params.csf.vmax, params.csf.beta, params.csf.trunc);
semilogx(f1D, csf1D, 'k-', 'LineWidth', 2.2); hold on;
xline(optSF, 'r--', sprintf('Optimal\n%.2f cpd', optSF), 'LineWidth', 1.4, 'LabelVerticalAlignment', 'bottom');
set(gca, 'XTick', [0.1 0.5 1 3 8 16], 'XTickLabel', {'0.1','0.5','1','3','8','16'});
xlim([0.05 16]); ylim([0 1.1]); grid on; axis square;
xlabel('Spatial Frequency (cpd)'); ylabel('CSF Gain');
title(sprintf('2a. Human Photopic CSF (f_{max}=%.1f cpd)', params.csf.fmax), 'FontSize', 11);

subplot(2, 4, 3);
imagesc(im6_csf, [-csfLim, csfLim]); axis image off; colormap(gca, gray(256));
title('2b. CSF-Filtered Display (N = 6)', 'FontSize', 11);

subplot(2, 4, 4);
imagesc(fZoom, fZoom, log10(max(P6_csf(zMask, zMask), 1e-6))); hold on;
plot(cx, cy, 'r--', 'LineWidth', 1.4);
axis xy image; colorbar;
xlabel('f_x (cpd)'); ylabel('f_y (cpd)');
title('3a. 2D CSF Power log_{10}(|FFT|^2\cdotCSF^2), N=6', 'FontSize', 11);

% Bottom row: N = 29
subplot(2, 4, 5);
imshow(im29, [0 1]);
title('1b. Test Display (N = 29)', 'FontSize', 11);

subplot(2, 4, 6);
imagesc(fZoom, fZoom, csf2D(zMask, zMask)); hold on;
plot(cx, cy, 'r--', 'LineWidth', 1.4);
axis xy image; colorbar;
xlabel('f_x (cpd)'); ylabel('f_y (cpd)');
title('2c. 2D CSF Filter (DC = 0)', 'FontSize', 11);

subplot(2, 4, 7);
imagesc(im29_csf, [-csfLim, csfLim]); axis image off; colormap(gca, gray(256));
title('2d. CSF-Filtered Display (N = 29)', 'FontSize', 11);

subplot(2, 4, 8);
imagesc(fZoom, fZoom, log10(max(P29_csf(zMask, zMask), 1e-6))); hold on;
plot(cx, cy, 'r--', 'LineWidth', 1.4);
axis xy image; colorbar;
xlabel('f_x (cpd)'); ylabel('f_y (cpd)');
title('3b. 2D CSF Power log_{10}(|FFT|^2\cdotCSF^2), N=29', 'FontSize', 11);

save_fig_png(fig2, fullfile(suppDir, 'supp_ideal_observer_steps1to3.png'));
close(fig2);

%% ========================================================================
%  FIGURE 3: Detailed Steps 4-5 (Spectral Differences -> MEG RSA)
%  ========================================================================
fig3 = figure('Visible', 'off', 'Name', 'Ideal Observer Steps 4-5 Detail', 'Color', 'w', 'Position', [70 110 1480 460]);

% Panel A: 2D Difference Spectrum (N=29 vs N=6, CSF-weighted)
subplot(1, 3, 1);
diff2D = abs(P29_csf(zMask, zMask) - P6_csf(zMask, zMask));
imagesc(fZoom, fZoom, log10(max(diff2D, 1e-6))); hold on;
plot(cx, cy, 'r--', 'LineWidth', 1.6);
axis xy image; colorbar;
xlabel('f_x (cpd)'); ylabel('f_y (cpd)');
title({'4a. 2D Spectral Difference', '|P_{CSF}(N=29) - P_{CSF}(N=6)|'}, 'FontSize', 11);

% Panel B: 1D Pairwise Spectral Differences & Peak Alignment
subplot(1, 3, 2); hold on;
for p = 1:nPairs
    plot(fSearch, opt.pairCurves(p, sfMask), 'Color', pCols(p, :), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s (%.2f cpd)', opt.pairLabels{p}, opt.pairPeaks_cpd(p)));
end
xline(optSF, 'r--', sprintf('Target = %.2f cpd', optSF), 'LineWidth', 2.0, ...
    'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.5 1 2 4 8 16], ...
    'XTickLabel', {'0.5','1','2','4','8','16'});
xlim([params.csf.minSF_cpd, 16]);
xlabel('Spatial Frequency (cpd, > 0.5 cpd)'); ylabel('|\Delta Power|');
title('4b. Pairwise Differences & Target SF', 'FontSize', 11);
legend('Location', 'southwest', 'FontSize', 8); grid on; axis square;

% Panel C: Step 5 - MEG RSA SF Tuning Curve (or 2.56 cpd RDM if MEG summary absent)
subplot(1, 3, 3);
if hasMegSum
    nB = numel(sfFreqs_cpd);
    eVals = sfTbl.Early_80_140ms_dR2(1:nB);
    lVals = sfTbl.Late_220_355ms_dR2(1:nB);
    eNum  = sfTbl.Early_80_140ms_dR2(strcmp(sfTbl.Band, 'Number_Ref'));
    lNum  = sfTbl.Late_220_355ms_dR2(strcmp(sfTbl.Band, 'Number_Ref'));
    hold on;
    plot(sfFreqs_cpd, eVals, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 2.0, ...
        'MarkerFaceColor', [0.00 0.45 0.74], 'DisplayName', 'Early (80-140 ms)');
    plot(sfFreqs_cpd, lVals, '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.0, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Late (220-355 ms)');
    yline(eNum, '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.2, 'DisplayName', 'Number (Early)');
    yline(lNum, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, 'DisplayName', 'Number (Late)');
    xline(optSF, 'r:', 'LineWidth', 1.5, 'DisplayName', sprintf('Ideal Obs (%.2f cpd)', optSF));
    set(gca, 'XScale', 'log');
    xlim([sfFreqs_cpd(1)*0.8, sfFreqs_cpd(end)*1.25]);
    xticks(sfFreqs_cpd); xticklabels(arrayfun(@(f) sprintf('%.2f', f), sfFreqs_cpd, 'UniformOutput', false));
    xlabel('Spatial Frequency (cpd)'); ylabel('Incremental \DeltaR^2 (above 4 non-num)');
    title('5. MEG RSA Tuning Peaks at Target SF', 'FontSize', 11);
    legend('Location', 'northwest', 'FontSize', 8); grid on; axis square;
else
    imagesc(rdm.SF_Bands{optBandIdx}); axis square; colorbar;
    title(sprintf('5. Target %.2f cpd RDM for RSA', sfFreqs_cpd(optBandIdx)), 'FontSize', 11);
    xlabel('Condition'); ylabel('Condition');
end

save_fig_png(fig3, fullfile(suppDir, 'supp_ideal_observer_steps4to5.png'));
close(fig3);
fprintf('Saved 3 ideal-observer walkthrough figures to:\n  %s\n', suppDir);

function save_fig_png(figHandle, filepath)
if exist('exportgraphics', 'file')
    exportgraphics(figHandle, filepath, 'Resolution', 150);
else
    saveas(figHandle, filepath);
end
end
