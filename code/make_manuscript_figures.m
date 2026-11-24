%% MAKE_MANUSCRIPT_FIGURES
% Generates all 7 publication-ready standalone 300-DPI panels for the
% revised 2-figure Matters Arising manuscript (N = 30), using Karami et al.'s
% exact feature terminology (Number, Average Item Area, Total Field Area,
% Total Surface Area, Density).
%
% Outputs (saved to FRF_model/results/panels_main/):
%   Figure 1 panels:
%     - fig1a_karami_replication.png    (Karami 5-predictor semipartial RSA, N=30)
%     - fig1b_stimuli_and_csf.png       (2x2 real stimulus displays: raw & CSF-filtered)
%     - fig1c_spectral_differences.png  (Clean pairwise CSF spectral differences, no title/cutoff)
%     - fig1d_ce_rdm.png                (32x32 RDM for critical 2.56 cpd contrast energy)
%   Figure 2 panels:
%     - fig2a_semipartial_ce_vs_number.png (Semipartial r: Number vs. Contrast Energy)
%     - fig2b_variance_partitioning.png    (Non-negative variance partitioning, clean legend)
%     - fig2c_karami_with_ce.png           (6-predictor RSA: Number & CE full contrast + SEM;
%                                           4 non-numeric covariates faded without error bars)

clear; close all; clc;

baseDir   = fileparts(mfilename('fullpath'));
resDir    = fullfile(baseDir, 'results');
mainDir   = fullfile(resDir, 'panels_main');
megRdmDir = fullfile(baseDir, 'meg_rdms');
if ~exist(mainDir, 'dir'), mkdir(mainDir); end

params = frf_default_params();
addpath(params.filterDir);

%% 1. Load ideal-observer, model RDMs, and all 30 subject MEG RDMs
load(fullfile(resDir, 'optimal_sf.mat'), 'opt');
M = load(fullfile(resDir, 'frf_rdms.mat'), 'rdm', 'sfFreqs_cpd', 'optBandSF_cpd');
rdm           = M.rdm;
optBandSF_cpd = M.optBandSF_cpd;
optSF         = opt.optimalSF_cpd;

utMask = triu(true(32), 1);
vNum   = zscore(rdm.Number(utMask));
vAIA   = zscore(rdm.ItemArea(utMask));
vTFA   = zscore(rdm.FieldArea(utMask));
vTSA   = zscore(rdm.SurfaceArea(utMask));
vDen   = zscore(rdm.Density(utMask));
vCE    = zscore(rdm.ContrastEnergy(utMask)); % Single 2.56 cpd band

% Exact Karami et al. (2026) predictor names
karamiLabels = {'Number', 'Total Field Area', 'Total Surface Area', 'Average Item Area', 'Density'};
X_karami5    = [vNum, vTFA, vTSA, vAIA, vDen];
X_nonNum     = [vAIA, vTFA, vTSA, vDen];
X_karami6    = [vNum, vCE, vTFA, vTSA, vAIA, vDen]; % 6-predictor model adding CE

% Precompute residuals for 5-predictor Karami RSA
res_k5 = zeros(size(X_karami5));
for j = 1:5
    res_k5(:, j) = residualize(X_karami5(:, j), X_karami5(:, setdiff(1:5, j)));
end

% Precompute residuals for 6-predictor RSA (Karami + 2.56 cpd CE)
res_k6 = zeros(size(X_karami6));
for j = 1:6
    res_k6(:, j) = residualize(X_karami6(:, j), X_karami6(:, setdiff(1:6, j)));
end

% Residuals for CE | 4 non-numeric covariates & 2-predictor NNLS VP
res_Num_nonNum = res_k5(:, 1);
res_CE_nonNum  = residualize(vCE, X_nonNum);
X_pair_r       = [res_Num_nonNum, res_CE_nonNum];

subjFiles = dir(fullfile(megRdmDir, 'sub-*_rdm_grad.mat'));
nSubj     = numel(subjFiles);
assert(nSubj > 0, 'No subject MEG RDM files found in %s.', megRdmDir);
fprintf('Loaded %d subject MEG RDM files.\n', nSubj);

tmp      = load(fullfile(megRdmDir, subjFiles(1).name));
times_ms = double(tmp.times_ms(:)');
nTimes   = numel(times_ms);

sr_k5        = zeros(nSubj, nTimes, 5);
sr_k6        = zeros(nSubj, nTimes, 6);
sr_CE_nonNum = zeros(nSubj, nTimes);
vp_shared    = zeros(nSubj, nTimes);
vp_uniqNum   = zeros(nSubj, nTimes);
vp_uniqCE    = zeros(nSubj, nTimes);

for s = 1:nSubj
    S = load(fullfile(megRdmDir, subjFiles(s).name));
    megRDM = double(S.rdm_grad);
    for t = 1:nTimes
        y = megRDM(:, :, t);
        y = y(utMask);
        if std(y) == 0, continue; end
        
        for j = 1:5
            sr_k5(s, t, j) = corr(y, res_k5(:, j));
        end
        for j = 1:6
            sr_k6(s, t, j) = corr(y, res_k6(:, j));
        end
        sr_Num = sr_k5(s, t, 1);
        sr_CE  = corr(y, res_CE_nonNum);
        sr_CE_nonNum(s, t) = sr_CE;
        
        % Non-negative 2-predictor variance partitioning
        dNum  = max(0, sr_Num)^2;
        dCE   = max(0, sr_CE)^2;
        y_r   = residualize(y, X_nonNum);
        ssTot = sum((y - mean(y)).^2);
        b_ols = X_pair_r \ y_r;
        if b_ols(1) >= 0 && b_ols(2) >= 0
            dBoth = max(0, (sum(y_r.^2) - sum((y_r - X_pair_r * b_ols).^2)) / ssTot);
        else
            dBoth = max(dNum, dCE);
        end
        vp_uniqNum(s, t) = max(0, dBoth - dCE);
        vp_uniqCE(s, t)  = max(0, dBoth - dNum);
        vp_shared(s, t)  = max(0, dNum + dCE - dBoth);
    end
end

tMask = times_ms >= -50 & times_ms <= 600;
tPlot = times_ms(tMask);

% Consistent color palette
colNum  = [0.00 0.00 0.00]; % Black
colCE   = [0.00 0.45 0.74]; % Blue
colTFA  = [0.85 0.33 0.10]; % Vermilion/Orange
colTSA  = [0.15 0.62 0.65]; % Teal
colAIA  = [0.47 0.67 0.19]; % Green
colDen  = [0.58 0.25 0.62]; % Purple

k5Cols  = [colNum; colTFA; colTSA; colAIA; colDen];
k5Lines = {'-', '-', '-', '-', '--'};

%% ========================================================================
%  FIGURE 1 PANELS (4 panels)
%  ========================================================================

% --- Panel 1a: Karami et al. 5-predictor replication ---
f1a = figure('Color', 'w', 'Position', [60 400 580 440]);
hold on;
shade_window([80 140],  [-0.10 0.34], [0.92 0.95 1.0]);
shade_window([220 355], [-0.10 0.34], [1.0 0.95 0.92]);
h1a = gobjects(1, 5);
for j = 1:5
    h1a(j) = plot_mean_sem(tPlot, sr_k5(:, tMask, j), k5Cols(j, :), k5Lines{j}, 2.0, karamiLabels{j});
end
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.08 0.33]);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Semipartial Pearson r', 'FontSize', 11);
title(sprintf('Karami et al. Replication (n = %d)', nSubj), 'FontSize', 11.5);
legend(h1a, 'Location', 'northeast', 'FontSize', 8.5, 'Box', 'on');
set(gca, 'FontSize', 10); grid on;
save_fig_png(f1a, fullfile(mainDir, 'fig1a_karami_replication.png'), 300);

% --- Panel 1b: 2x2 Real Stimulus Displays (Raw top row, CSF-filtered bottom row) ---
% Matched Total Surface Area mid-range pair from Karami et al.'s 32-condition grid:
%   N = 10 (Average Item Area = 0.20 vd^2, TSA = 2.00 vd^2, Total Field Area = 13.5 vd)
%   N = 17 (Average Item Area = 0.12 vd^2, TSA = 2.04 vd^2, Total Field Area = 13.5 vd)
rng(params.stim.rngSeed);
exTFA = 13.5;
im10  = generate_karami_stimulus(10, 0.20, exTFA, params.stim);
im17  = generate_karami_stimulus(17, 0.12, exTFA, params.stim);

imSize = params.stim.imSize_px;
imDeg  = params.stim.imSize_deg;
[~, radDist] = polarDistFun(imSize);
radCPD = radDist / imDeg;
csf2D  = createCSF(max(radCPD, 1e-4), params.csf.fmax, params.csf.vmax, params.csf.beta, params.csf.trunc);
csf2D(radCPD == 0) = 0;

im10_csf = real(ifft2(ifftshift(fftshift(fft2(im10 - mean(im10(:)))) .* csf2D)));
im17_csf = real(ifft2(ifftshift(fftshift(fft2(im17 - mean(im17(:)))) .* csf2D)));
csfLim   = max([abs(im10_csf(:)); abs(im17_csf(:))]) * 0.85;

f1b = figure('Color', 'w', 'Position', [660 400 500 460]);
tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(tlo);
imshow(im10, [0 1]);
title('N = 10 (0.20 vd^2, 13.5\circ)', 'FontSize', 9.5);

nexttile(tlo);
imshow(im17, [0 1]);
title('N = 17 (0.12 vd^2, 13.5\circ)', 'FontSize', 9.5);

nexttile(tlo);
imagesc(im10_csf, [-csfLim, csfLim]); axis image off; colormap(gca, gray(256));
title('CSF-Filtered (N = 10)', 'FontSize', 9.5);

nexttile(tlo);
imagesc(im17_csf, [-csfLim, csfLim]); axis image off; colormap(gca, gray(256));
title('CSF-Filtered (N = 17)', 'FontSize', 9.5);

save_fig_png(f1b, fullfile(mainDir, 'fig1b_stimuli_and_csf.png'), 300);

% --- Panel 1c: Clean Pairwise Spectral Differences (No title, no 0.5 cpd cutoff line) ---
% Pick the 3 adjacent numerosity comparisons from opt.pairLabels:
%   1: N=6 vs 10, 4: N=10 vs 17, 6: N=17 vs 29
pickPairs = [1, 4, 6];
pairCols  = [0.18 0.35 0.65; 0.00 0.55 0.55; 0.85 0.45 0.10];

sfMask  = opt.fBins_cpd >= 0.5;
fSearch = opt.fBins_cpd(sfMask);

f1c = figure('Color', 'w', 'Position', [1180 400 520 440]);
hold on;
h1c = gobjects(1, numel(pickPairs) + 1);
for idx = 1:numel(pickPairs)
    p = pickPairs(idx);
    h1c(idx) = plot(fSearch, opt.pairCurves(p, sfMask), 'Color', pairCols(idx, :), ...
        'LineWidth', 2.2, 'DisplayName', opt.pairLabels{p});
end
h1c(end) = xline(optSF, 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Peak (%.2f cpd)', optSF));
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XTick', [0.5 1 2 4 8 16], 'XTickLabel', {'0.5','1','2','4','8','16'}, 'FontSize', 10);
xlim([0.5, 16]);
xlabel('Spatial Frequency (cpd)', 'FontSize', 11);
ylabel('|\Delta CSF-Weighted Power|', 'FontSize', 11);
legend(h1c, 'Location', 'southwest', 'FontSize', 9, 'Box', 'on');
grid on;
save_fig_png(f1c, fullfile(mainDir, 'fig1c_spectral_differences.png'), 300);

% --- Panel 1d: 32x32 RDM for Critical 2.56 cpd Contrast Energy ---
f1d = figure('Color', 'w', 'Position', [100 40 480 440]);
imagesc(rdm.ContrastEnergy);
axis square;
cb = colorbar;
cb.Label.String = 'Distance (log)';
cb.Label.FontSize = 10;
set(gca, 'XTick', [1 8 16 24 32], 'YTick', [1 8 16 24 32], 'FontSize', 10);
xlabel('Stimulus Condition', 'FontSize', 11);
ylabel('Stimulus Condition', 'FontSize', 11);
title(sprintf('Contrast Energy RDM (%.2f cpd)', optBandSF_cpd), 'FontSize', 11.5);
save_fig_png(f1d, fullfile(mainDir, 'fig1d_ce_rdm.png'), 300);

%% ========================================================================
%  FIGURE 2 PANELS (3 panels)
%  ========================================================================

% --- Panel 2a: Semipartial r — Number vs Contrast Energy (clean legend) ---
f2a = figure('Color', 'w', 'Position', [600 40 580 440]);
hold on;
shade_window([80 140],  [-0.06 0.22], [0.92 0.95 1.0]);
shade_window([220 355], [-0.06 0.22], [1.0 0.95 0.92]);
h2a = gobjects(1, 2);
h2a(1) = plot_mean_sem(tPlot, sr_k5(:, tMask, 1),     colNum, '-', 2.2, 'Number');
h2a(2) = plot_mean_sem(tPlot, sr_CE_nonNum(:, tMask), colCE,  '-', 2.2, sprintf('Contrast Energy (%.2f cpd)', optBandSF_cpd));
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.05 0.21]);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Semipartial Pearson r', 'FontSize', 11);
title('Number vs. Contrast Energy', 'FontSize', 11.5);
legend(h2a, 'Location', 'northeast', 'FontSize', 9.5, 'Box', 'on');
set(gca, 'FontSize', 10); grid on;
save_fig_png(f2a, fullfile(mainDir, 'fig2a_semipartial_ce_vs_number.png'), 300);

% --- Panel 2b: Non-Negative Variance Partitioning (clean legend) ---
f2b = figure('Color', 'w', 'Position', [1200 40 580 440]);
hold on;
shade_window([80 140],  [-0.002 0.045], [0.92 0.95 1.0]);
shade_window([220 355], [-0.002 0.045], [1.0 0.95 0.92]);
muShrd = mean(vp_shared(:, tMask), 1);
hSh = area(tPlot, max(0, muShrd), 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.50, ...
    'EdgeColor', [0.48 0.48 0.48], 'LineWidth', 1.2, 'DisplayName', 'Shared variance');
hUN = plot_mean_sem(tPlot, vp_uniqNum(:, tMask), colNum, '-', 2.2, 'Unique Number');
hUE = plot_mean_sem(tPlot, vp_uniqCE(:, tMask),  colCE,  '-', 2.2, 'Unique Contrast Energy');
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.002 0.042]);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Variance Explained (\DeltaR^2)', 'FontSize', 11);
title('Variance Partitioning', 'FontSize', 11.5);
legend([hSh, hUN, hUE], 'Location', 'northeast', 'FontSize', 9.5, 'Box', 'on');
set(gca, 'FontSize', 10); grid on;
save_fig_png(f2b, fullfile(mainDir, 'fig2b_variance_partitioning.png'), 300);

% --- Panel 2c: Karami et al. RSA with Contrast Energy Included (6 predictors) ---
% Non-numeric predictors faded into background without error bars;
% Number and Contrast Energy (2.56 cpd) at full contrast with SEM shading.
f2c = figure('Color', 'w', 'Position', [300 200 580 440]);
hold on;
shade_window([80 140],  [-0.10 0.34], [0.92 0.95 1.0]);
shade_window([220 355], [-0.10 0.34], [1.0 0.95 0.92]);

% Plot 4 non-numeric covariates first (faded, no error bars)
fadeFactor = 0.72; % blend 72% toward white so they recede cleanly
bgIdx   = [3, 4, 5, 6]; % TFA, TSA, AIA, Density in X_karami6
bgCols  = [colTFA; colTSA; colAIA; colDen];
bgLines = {'-', '-', '-', '--'};
bgNames = {'Total Field Area', 'Total Surface Area', 'Average Item Area', 'Density'};
hBg = gobjects(1, 4);
for k = 1:4
    fadedCol = bgCols(k, :) * (1 - fadeFactor) + [1 1 1] * fadeFactor;
    muBg = mean(sr_k6(:, tMask, bgIdx(k)), 1);
    hBg(k) = plot(tPlot, muBg, 'Color', fadedCol, 'LineStyle', bgLines{k}, ...
        'LineWidth', 1.4, 'DisplayName', bgNames{k});
end

% Plot Number and Contrast Energy at full contrast with SEM bands
hNum6 = plot_mean_sem(tPlot, sr_k6(:, tMask, 1), colNum, '-', 2.3, 'Number');
hCE6  = plot_mean_sem(tPlot, sr_k6(:, tMask, 2), colCE,  '-', 2.3, sprintf('Contrast Energy (%.2f cpd)', optBandSF_cpd));

xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.08 0.33]);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Semipartial Pearson r', 'FontSize', 11);
title('Full Model Including Contrast Energy', 'FontSize', 11.5);
legend([hNum6, hCE6, hBg], 'Location', 'northeast', 'FontSize', 8.5, 'Box', 'on');
set(gca, 'FontSize', 10); grid on;
save_fig_png(f2c, fullfile(mainDir, 'fig2c_karami_with_ce.png'), 300);

fprintf('\nSaved all 7 standalone 300-DPI manuscript panels to:\n  %s\n', mainDir);

%% Helper functions
function r = residualize(x, Z)
Zd = [ones(size(Z, 1), 1), Z];
r  = x - Zd * (Zd \ x);
end

function h = plot_mean_sem(t, dat, col, lineStyle, lw, dispName)
mu = mean(dat, 1);
if size(dat, 1) >= 2
    sem = std(dat, 0, 1) / sqrt(size(dat, 1));
    fill([t, fliplr(t)], [mu + sem, fliplr(mu - sem)], col, ...
        'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
h = plot(t, mu, 'Color', col, 'LineStyle', lineStyle, 'LineWidth', lw, 'DisplayName', dispName);
end

function shade_window(xRange, yRange, col)
patch([xRange(1) xRange(2) xRange(2) xRange(1)], ...
      [yRange(1) yRange(1) yRange(2) yRange(2)], ...
      col, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

function save_fig_png(figHandle, filepath, dpi)
if nargin < 3, dpi = 300; end
if exist('exportgraphics', 'file')
    exportgraphics(figHandle, filepath, 'Resolution', dpi);
else
    saveas(figHandle, filepath);
end
end
