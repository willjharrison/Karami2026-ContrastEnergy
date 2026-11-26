%% RUN_MEG_VP_TEST
% Stage B MEG analysis testing whether local contrast energy at the single
% CSF-limited ideal-observer spatial frequency accounts for Karami et al.'s
% (2026) MEG Number effect, with supplementary multi-band SF tuning curves.
%
% Uses non-negative (positive-fit) incremental R^2 and signed semipartial r
% so negative suppressor weights are never credited as positive model fit.
%
% Outputs:
%   results/panels_main/ (single ideal-observer SF; displayed on screen, 300 DPI)
%     - main_1_semipartial_ce_vs_number.png
%     - main_2_variance_partitioning.png
%     - main_3_matched_5pred_diff.png
%   results/panels_supplement/ (saved silently, 300 DPI)
%     - supp_1_karami_replication.png
%     - supp_2_zero_order_fits.png
%     - supp_3_sf_band_timecourses.png
%     - supp_4_sf_time_heatmap.png
%     - supp_5_sf_tuning_early_vs_late.png
%   results/
%     - meg_vp_summary.csv
%     - sf_band_summary.csv

clear; close all; clc;

baseDir   = fileparts(mfilename('fullpath'));
resDir    = fullfile(baseDir, 'results');
mainDir   = fullfile(resDir, 'panels_main');
suppDir   = fullfile(resDir, 'panels_supplement');
megRdmDir = fullfile(baseDir, 'meg_rdms');

if ~exist(mainDir, 'dir'), mkdir(mainDir); end
if ~exist(suppDir, 'dir'), mkdir(suppDir); end

%% 1. Load Contrast Energy & Karami model RDMs
modelFile = fullfile(resDir, 'frf_rdms.mat');
assert(exist(modelFile, 'file') == 2, 'Run run_frf_rdm_simulation.m first to create frf_rdms.mat');
M = load(modelFile, 'rdm', 'sfFreqs_cpd', 'sfBandNames', 'optBandIdx', 'optBandSF_cpd');
rdm         = M.rdm;
sfFreqs_cpd = M.sfFreqs_cpd;
nBands      = numel(sfFreqs_cpd);
if isfield(M, 'optBandSF_cpd')
    optBandSF_cpd = M.optBandSF_cpd;
else
    optBandSF_cpd = sfFreqs_cpd(ceil(nBands / 2));
end

utMask    = triu(true(32), 1);
vNum      = zscore(rdm.Number(utMask));
vAIA      = zscore(rdm.ItemArea(utMask));
vTFA      = zscore(rdm.FieldArea(utMask));
vTSA      = zscore(rdm.SurfaceArea(utMask));
vDen      = zscore(rdm.Density(utMask));
vCE       = zscore(rdm.ContrastEnergy(utMask));           % Single ideal-observer SF band
vCE_multi = zscore(rdm.ContrastEnergy_MultiBand(utMask)); % Supplementary multi-band bank
vCont     = zscore(rdm.Contour(utMask));
vLum      = zscore(rdm.MeanLuminance(utMask));

karamiNames = {'Number', 'Total Field Area', 'Total Surface Area', 'Average Item Area', 'Density'};
X_karami    = [vNum, vTFA, vTSA, vAIA, vDen];
X_nonNum    = [vAIA, vTFA, vTSA, vDen];

V_bands = zeros(sum(utMask(:)), nBands);
for b = 1:nBands
    V_bands(:, b) = zscore(rdm.SF_Bands{b}(utMask));
end

%% 2. Load all available subject MEG RDMs
subjFiles = dir(fullfile(megRdmDir, 'sub-*_rdm_grad.mat'));
assert(~isempty(subjFiles), 'No sub-*_rdm_grad.mat files found in %s.', megRdmDir);
nSubj = numel(subjFiles);
fprintf('Loaded %d subject MEG RDM file(s) from %s\n', nSubj, megRdmDir);

tmp = load(fullfile(megRdmDir, subjFiles(1).name));
times_ms = double(tmp.times_ms(:)');
nTimes   = numel(times_ms);

% Precompute residuals (orthogonalized to 4 non-numeric covariates)
res_karami = zeros(size(X_karami));
for j = 1:size(X_karami, 2)
    otherCols = setdiff(1:size(X_karami, 2), j);
    res_karami(:, j) = residualize(X_karami(:, j), X_karami(:, otherCols));
end
res_Num_nonNum     = res_karami(:, 1);
res_CE_nonNum      = residualize(vCE, X_nonNum);
res_CEmulti_nonNum = residualize(vCE_multi, X_nonNum);
res_Cont_nonNum    = residualize(vCont, X_nonNum);
res_Lum_nonNum     = residualize(vLum, X_nonNum);

% Unique residuals (controlling for 4 non-numeric covariates + competing model)
res_Num_afterCE = residualize(vNum, [X_nonNum, vCE]);
res_CE_afterNum = residualize(vCE,  [X_nonNum, vNum]);

res_bands = zeros(size(V_bands));
for b = 1:nBands
    res_bands(:, b) = residualize(V_bands(:, b), X_nonNum);
end

% Design matrix for exact 2-predictor NNLS variance partitioning
X_pair_r = [res_Num_nonNum, res_CE_nonNum];

% Preallocate subject x time arrays
sr_karami    = zeros(nSubj, nTimes, numel(karamiNames));
sr_CE_nonNum = zeros(nSubj, nTimes);
sr_uniqNum   = zeros(nSubj, nTimes); % Signed: Number | 4 non-num + CE
sr_uniqCE    = zeros(nSubj, nTimes); % Signed: CE | 4 non-num + Number

r0_Num   = zeros(nSubj, nTimes);
r0_CE    = zeros(nSubj, nTimes);
r0_Cont  = zeros(nSubj, nTimes);
r0_Lum   = zeros(nSubj, nTimes);

r2diff_CE   = zeros(nSubj, nTimes); % Positive-fit dR^2(CE) - dR^2(Num)
r2diff_Cont = zeros(nSubj, nTimes); % Positive-fit dR^2(Contour) - dR^2(Num)
r2diff_Lum  = zeros(nSubj, nTimes); % Positive-fit dR^2(MeanLuminance) - dR^2(Num)

vp_shared  = zeros(nSubj, nTimes);
vp_uniqNum = zeros(nSubj, nTimes);
vp_uniqCE  = zeros(nSubj, nTimes);

dR2_Num      = zeros(nSubj, nTimes);
dR2_CE       = zeros(nSubj, nTimes);
dR2_CE_multi = zeros(nSubj, nTimes);
dR2_bands    = zeros(nSubj, nTimes, nBands);

%% 3. Run time-resolved regressions per subject
for s = 1:nSubj
    S = load(fullfile(megRdmDir, subjFiles(s).name));
    megRDM = double(S.rdm_grad);
    
    for t = 1:nTimes
        y = megRDM(:, :, t);
        y = y(utMask);
        if std(y) == 0, continue; end
        
        % Signed semipartial correlations (exact Karami et al. metric)
        for j = 1:numel(karamiNames)
            sr_karami(s, t, j) = corr(y, res_karami(:, j));
        end
        sr_Num  = sr_karami(s, t, 1);
        sr_CE   = corr(y, res_CE_nonNum);
        sr_CEm  = corr(y, res_CEmulti_nonNum);
        sr_Cont = corr(y, res_Cont_nonNum);
        sr_Lum  = corr(y, res_Lum_nonNum);
        
        sr_CE_nonNum(s, t) = sr_CE;
        sr_uniqNum(s, t)   = corr(y, res_Num_afterCE);
        sr_uniqCE(s, t)    = corr(y, res_CE_afterNum);
        
        % Zero-order correlations
        r0_Num(s, t)  = corr(y, vNum);
        r0_CE(s, t)   = corr(y, vCE);
        r0_Cont(s, t) = corr(y, vCont);
        r0_Lum(s, t)  = corr(y, vLum);
        
        % Positive-fit incremental R^2 above 4 non-numeric covariates:
        % For a single predictor v, positive-fit dR^2 = max(0, semipartial r)^2.
        dNum  = max(0, sr_Num)^2;
        dCE   = max(0, sr_CE)^2;
        dCEm  = max(0, sr_CEm)^2;
        dCont = max(0, sr_Cont)^2;
        dLum  = max(0, sr_Lum)^2;
        
        dR2_Num(s, t)      = dNum;
        dR2_CE(s, t)       = dCE;
        dR2_CE_multi(s, t) = dCEm;
        
        r2diff_CE(s, t)   = dCE   - dNum;
        r2diff_Cont(s, t) = dCont - dNum;
        r2diff_Lum(s, t)  = dLum  - dNum;
        
        % Exact 2-predictor non-negative variance partitioning (eliminates
        % negative-beta suppressor artifacts between collinear predictors)
        y_r    = residualize(y, X_nonNum);
        ssTot  = sum((y - mean(y)).^2);
        b_ols  = X_pair_r \ y_r;
        if b_ols(1) >= 0 && b_ols(2) >= 0
            ssRes_0    = sum(y_r.^2);
            ssRes_both = sum((y_r - X_pair_r * b_ols).^2);
            dBoth      = max(0, (ssRes_0 - ssRes_both) / ssTot);
        else
            % Convex active-set boundary: NNLS optimum sets the negative beta to 0
            dBoth = max(dNum, dCE);
        end
        
        vp_uniqNum(s, t) = max(0, dBoth - dCE);
        vp_uniqCE(s, t)  = max(0, dBoth - dNum);
        vp_shared(s, t)  = max(0, dNum + dCE - dBoth);
        
        % Supplementary individual SF bands (positive-fit dR^2)
        for b = 1:nBands
            sr_b = corr(y, res_bands(:, b));
            dR2_bands(s, t, b) = max(0, sr_b)^2;
        end
    end
end

%% 4. Print & save summaries
earlyMask = times_ms >= 115 & times_ms <= 140;
lateMask  = times_ms >= 220 & times_ms <= 355;

fprintf('\n=== Primary Analysis: Single Ideal-Observer SF (%.2f cpd) vs. Number (N = %d) ===\n', optBandSF_cpd, nSubj);
eDiff_s  = mean(r2diff_CE(:, earlyMask), 2);
lDiff_s  = mean(r2diff_CE(:, lateMask), 2);
eUnNum_s = mean(sr_uniqNum(:, earlyMask), 2);
lUnNum_s = mean(sr_uniqNum(:, lateMask), 2);
eUnCE_s  = mean(sr_uniqCE(:, earlyMask), 2);
lUnCE_s  = mean(sr_uniqCE(:, lateMask), 2);

if nSubj >= 3
    [~, pE, ~, stE]     = ttest(eDiff_s);
    [~, pL, ~, stL]     = ttest(lDiff_s);
    [~, pUN_E, ~, stUNE] = ttest(eUnNum_s);
    [~, pUN_L, ~, stUNL] = ttest(lUnNum_s);
    [~, pUC_E, ~, stUCE] = ttest(eUnCE_s);
    [~, pUC_L, ~, stUCL] = ttest(lUnCE_s);
    statE = sprintf(' | t(%d)=%+.2f, p=%.4f', stE.df, stE.tstat, pE);
    statL = sprintf(' | t(%d)=%+.2f, p=%.4f', stL.df, stL.tstat, pL);
else
    pE = NaN; pL = NaN; pUN_E = NaN; pUN_L = NaN; pUC_E = NaN; pUC_L = NaN;
    statE = ''; statL = '';
end

eSh  = mean(mean(vp_shared(:, earlyMask)));  lSh  = mean(mean(vp_shared(:, lateMask)));
eUN  = mean(mean(vp_uniqNum(:, earlyMask))); lUN  = mean(mean(vp_uniqNum(:, lateMask)));
eUC  = mean(mean(vp_uniqCE(:, earlyMask)));  lUC  = mean(mean(vp_uniqCE(:, lateMask)));
eNum = mean(mean(dR2_Num(:, earlyMask)));    lNum = mean(mean(dR2_Num(:, lateMask)));
eCE  = mean(mean(dR2_CE(:, earlyMask)));     lCE  = mean(mean(dR2_CE(:, lateMask)));
ePct = 100 * eSh / max(eNum, eps);           lPct = 100 * lSh / max(lNum, eps);

fprintf('Early (115-140 ms):  dR^2(CE - Num) = %+.4f%s | Shared=%.4f (%.1f%% of Num), UniqNum=%.4f, UniqCE=%.4f | CE/Num=%.1f%%\n', ...
    mean(eDiff_s), statE, eSh, ePct, eUN, eUC, 100 * eCE / max(eNum, eps));
fprintf('  Signed unique r:  Num|CE = %+.4f (t(%d)=%+.2f, p=%.4f) | CE|Num = %+.4f (t(%d)=%+.2f, p=%.4f)\n', ...
    mean(eUnNum_s), stUNE.df, stUNE.tstat, pUN_E, mean(eUnCE_s), stUCE.df, stUCE.tstat, pUC_E);
fprintf('Late  (220-355 ms): dR^2(CE - Num) = %+.4f%s | Shared=%.4f (%.1f%% of Num), UniqNum=%.4f, UniqCE=%.4f | CE/Num=%.1f%%\n', ...
    mean(lDiff_s), statL, lSh, lPct, lUN, lUC, 100 * lCE / max(lNum, eps));
fprintf('  Signed unique r:  Num|CE = %+.4f (t(%d)=%+.2f, p=%.4f) | CE|Num = %+.4f (t(%d)=%+.2f, p=%.4f)\n', ...
    mean(lUnNum_s), stUNL.df, stUNL.tstat, pUN_L, mean(lUnCE_s), stUCL.df, stUCL.tstat, pUC_L);
fprintf('MeanLuminance ctrl: Early dR^2(Lum - Num) = %+.4f | Late dR^2(Lum - Num) = %+.4f\n\n', ...
    mean(mean(r2diff_Lum(:, earlyMask))), mean(mean(r2diff_Lum(:, lateMask))));

vpTbl = table({'Early_115_140ms'; 'Late_220_355ms'}, ...
    [mean(eDiff_s); mean(lDiff_s)], [pE; pL], ...
    [eSh; lSh], [ePct; lPct], [eUN; lUN], [eUC; lUC], [eNum; lNum], [eCE; lCE], ...
    [mean(eUnNum_s); mean(lUnNum_s)], [pUN_E; pUN_L], [mean(eUnCE_s); mean(lUnCE_s)], [pUC_E; pUC_L], ...
    'VariableNames', {'Window', 'R2diff_CE_minus_Num', 'p_val', 'Shared_R2', 'Shared_Pct_of_Num', ...
    'UniqNum_R2', 'UniqCE_R2', 'Num_dR2', 'CE_dR2', 'SignedUniqNum_r', 'p_UniqNum', 'SignedUniqCE_r', 'p_UniqCE'});
writetable(vpTbl, fullfile(resDir, 'meg_vp_summary.csv'));

fprintf('=== Supplementary Multi-Band SF Incremental R^2 (above 4 non-num) ===\n');
sfRows = cell(nBands + 2, 3);
for b = 1:nBands
    eB = mean(mean(dR2_bands(:, earlyMask, b)));
    lB = mean(mean(dR2_bands(:, lateMask,  b)));
    fprintf('  Band %d (%6.2f cpd) | Early (115-140 ms): %.4f | Late (220-355 ms): %.4f\n', ...
        b, sfFreqs_cpd(b), eB, lB);
    sfRows(b, :) = {sprintf('%.2f cpd', sfFreqs_cpd(b)), eB, lB};
end
eAll = mean(mean(dR2_CE_multi(:, earlyMask))); lAll = mean(mean(dR2_CE_multi(:, lateMask)));
fprintf('  %-18s | Early (115-140 ms): %.4f | Late (220-355 ms): %.4f\n', 'Full CE Bank', eAll, lAll);
fprintf('  %-18s | Early (115-140 ms): %.4f | Late (220-355 ms): %.4f\n', 'Number (ref)', eNum, lNum);
sfRows(nBands + 1, :) = {'Full_CE_Bank', eAll, lAll};
sfRows(nBands + 2, :) = {'Number_Ref',   eNum, lNum};
writetable(cell2table(sfRows, 'VariableNames', {'Band', 'Early_115_140ms_dR2', 'Late_220_355ms_dR2'}), ...
    fullfile(resDir, 'sf_band_summary.csv'));

%% 5. MAIN PANELS (Single Ideal-Observer SF — displayed on screen, 300 DPI)
tMask = times_ms >= -50 & times_ms <= 600;
tPlot = times_ms(tMask);
panPos = [100 120 640 460];

% Main Panel 1: Semipartial r — Number vs Contrast Energy (single ideal-observer SF)
fM1 = figure('Name', 'Main 1: Semipartial r (Number vs CE)', 'Color', 'w', 'Position', panPos);
hold on;
shade_window([115 140],  [-0.06 0.22], [0.92 0.95 1.0]);
shade_window([220 355], [-0.06 0.22], [1.0 0.95 0.92]);
hM1 = gobjects(1, 2);
hM1(1) = plot_mean_sem(tPlot, sr_karami(:, tMask, 1), [0 0 0],          '-', 2.2, 'Number | 4 non-num');
hM1(2) = plot_mean_sem(tPlot, sr_CE_nonNum(:, tMask), [0.00 0.45 0.74], '-', 2.2, sprintf('Contrast Energy (%.2f cpd) | 4 non-num', optBandSF_cpd));
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.05 0.21]);
xlabel('Time (ms)'); ylabel('Semipartial Pearson r');
title(sprintf('Semipartial r: Number vs. Contrast Energy (N=%d)', nSubj));
legend(hM1, 'Location', 'northeast'); grid on;
save_fig_png(fM1, fullfile(mainDir, 'main_1_semipartial_ce_vs_number.png'), 300);

% Main Panel 2: Non-Negative Variance Partitioning (Shared vs Unique)
fM2 = figure('Name', 'Main 2: Variance Partitioning', 'Color', 'w', 'Position', panPos + [50 -40 0 0]);
hold on;
shade_window([115 140],  [-0.002 0.045], [0.92 0.95 1.0]);
shade_window([220 355], [-0.002 0.045], [1.0 0.95 0.92]);
muShrd = mean(vp_shared(:, tMask), 1);
hSh = area(tPlot, max(0, muShrd), 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.50, ...
    'EdgeColor', [0.50 0.50 0.50], 'LineWidth', 1.2, 'DisplayName', 'Shared (Number \cap CE)');
hUN = plot_mean_sem(tPlot, vp_uniqNum(:, tMask), [0 0 0],          '-', 2.2, 'Unique Number (after CE)');
hUE = plot_mean_sem(tPlot, vp_uniqCE(:, tMask),  [0.00 0.45 0.74], '-', 2.2, sprintf('Unique CE %.2f cpd (after Number)', optBandSF_cpd));
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.002 0.045]);
xlabel('Time (ms)'); ylabel('Variance Explained (\DeltaR^2)');
title('Non-Negative Variance Partitioning: Number vs. Contrast Energy');
legend([hSh, hUN, hUE], 'Location', 'northeast'); grid on;
save_fig_png(fM2, fullfile(mainDir, 'main_2_variance_partitioning.png'), 300);

% Main Panel 3: Matched 5-Predictor Comparison (Positive-fit dR^2_Model - dR^2_Number)
fM3 = figure('Name', 'Main 3: Matched 5-Predictor Comparison', 'Color', 'w', 'Position', panPos + [100 -80 0 0]);
hold on;
shade_window([115 140],  [-0.035 0.025], [0.92 0.95 1.0]);
shade_window([220 355], [-0.035 0.025], [1.0 0.95 0.92]);
hM3 = gobjects(1, 3);
hM3(1) = plot_mean_sem(tPlot, r2diff_CE(:, tMask),   [0.00 0.45 0.74], '-',  2.2, sprintf('Contrast Energy (%.2f cpd) - Number', optBandSF_cpd));
hM3(2) = plot_mean_sem(tPlot, r2diff_Cont(:, tMask), [0.20 0.65 0.20], ':',  1.9, 'Contour - Number');
hM3(3) = plot_mean_sem(tPlot, r2diff_Lum(:, tMask),  [0.55 0.55 0.55], '--', 1.6, 'Mean Luminance - Number');
xline(0, 'k:'); yline(0, 'k-', 'LineWidth', 1.0);
xlim([-50 600]); ylim([-0.035 0.020]);
xlabel('Time (ms)'); ylabel('\DeltaR^2_{Model} - \DeltaR^2_{Number}');
title('Matched 5-Predictor Comparison (>0 Model Wins)');
legend(hM3, 'Location', 'southeast'); grid on;
save_fig_png(fM3, fullfile(mainDir, 'main_3_matched_5pred_diff.png'), 300);

%% 6. SUPPLEMENTARY PANELS (Saved silently to results/panels_supplement/, 300 DPI)
% Supp 1: Karami 5-predictor replication
fS1 = figure('Visible', 'off', 'Color', 'w', 'Position', panPos);
hold on;
shade_window([115 140],  [-0.12 0.35], [0.92 0.95 1.0]);
shade_window([220 355], [-0.12 0.35], [1.0 0.95 0.92]);
kColors = [0 0 0; 0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56];
kLines  = {'-', '-', '-', '-', '--'};
hS1     = gobjects(1, numel(karamiNames));
for j = 1:numel(karamiNames)
    hS1(j) = plot_mean_sem(tPlot, sr_karami(:, tMask, j), kColors(j, :), kLines{j}, 2.0, karamiNames{j});
end
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.10 0.34]);
xlabel('Time (ms)'); ylabel('Semipartial Pearson r');
title(sprintf('Karami et al. Replication (N=%d)', nSubj));
legend(hS1, 'Location', 'northeast'); grid on;
save_fig_png(fS1, fullfile(suppDir, 'supp_1_karami_replication.png'), 300);
close(fS1);

% Supp 2: Zero-order fits
fS2 = figure('Visible', 'off', 'Color', 'w', 'Position', panPos);
hold on;
shade_window([115 140],  [-0.18 0.45], [0.92 0.95 1.0]);
shade_window([220 355], [-0.18 0.45], [1.0 0.95 0.92]);
hS2 = gobjects(1, 4);
hS2(1) = plot_mean_sem(tPlot, r0_Num(:, tMask),  [0 0 0],          '-',  2.2, 'Number');
hS2(2) = plot_mean_sem(tPlot, r0_CE(:, tMask),   [0.00 0.45 0.74], '-',  2.0, sprintf('Contrast Energy (%.2f cpd)', optBandSF_cpd));
hS2(3) = plot_mean_sem(tPlot, r0_Cont(:, tMask), [0.20 0.65 0.20], ':',  1.9, 'Contour');
hS2(4) = plot_mean_sem(tPlot, r0_Lum(:, tMask),  [0.55 0.55 0.55], '--', 1.6, 'Mean Luminance');
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([-0.18 0.44]);
xlabel('Time (ms)'); ylabel('Zero-Order Pearson r');
title(sprintf('Zero-Order Fit to MEG RDMs (N=%d)', nSubj));
legend(hS2, 'Location', 'northeast'); grid on;
save_fig_png(fS2, fullfile(suppDir, 'supp_2_zero_order_fits.png'), 300);
close(fS2);

% Supp 3: Individual SF band timecourses
sfCols = parula(nBands + 1);
fS3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 120 820 440]);
hold on;
shade_window([115 140],  [0 0.05], [0.92 0.95 1.0]);
shade_window([220 355], [0 0.05], [1.0 0.95 0.92]);
hS3 = gobjects(1, nBands + 1);
for b = 1:nBands
    hS3(b) = plot_mean_sem(tPlot, dR2_bands(:, tMask, b), sfCols(b, :), '-', 1.8, ...
        sprintf('%.2f cpd', sfFreqs_cpd(b)));
end
hS3(nBands + 1) = plot_mean_sem(tPlot, dR2_Num(:, tMask), [0 0 0], ':', 2.0, 'Number (ref)');
xline(0, 'k:'); yline(0, 'k:');
xlim([-50 600]); ylim([0 0.05]);
xlabel('Time (ms)'); ylabel('Incremental \DeltaR^2 (above 4 non-num)');
title(sprintf('Individual SF Bands Across Time (N=%d)', nSubj));
legend(hS3, 'Location', 'eastoutside'); grid on;
save_fig_png(fS3, fullfile(suppDir, 'supp_3_sf_band_timecourses.png'), 300);
close(fS3);

% Supp 4: SF x Time 2D heatmap
fS4 = figure('Visible', 'off', 'Color', 'w', 'Position', panPos);
sfTimeMap = squeeze(mean(dR2_bands(:, tMask, :), 1))';
imagesc(tPlot, 1:nBands, sfTimeMap);
set(gca, 'YDir', 'normal', 'YTick', 1:nBands, ...
    'YTickLabel', arrayfun(@(f) sprintf('%.2f', f), sfFreqs_cpd, 'UniformOutput', false));
xline(0, 'w:', 'LineWidth', 1.2);
xlabel('Time (ms)'); ylabel('Spatial Frequency (cpd)');
title('SF \times Time Map (\DeltaR^2 above 4 non-num)');
colorbar;
save_fig_png(fS4, fullfile(suppDir, 'supp_4_sf_time_heatmap.png'), 300);
close(fS4);

% Supp 5: Early vs Late SF tuning curve (showing peak at ideal-observer SF)
fS5 = figure('Visible', 'off', 'Color', 'w', 'Position', panPos);
hold on;
eBand_mu = squeeze(mean(mean(dR2_bands(:, earlyMask, :), 2), 1));
lBand_mu = squeeze(mean(mean(dR2_bands(:, lateMask,  :), 2), 1));
plot(sfFreqs_cpd, eBand_mu, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, 'MarkerFaceColor', [0.00 0.45 0.74], 'DisplayName', 'Early (115-140 ms)');
plot(sfFreqs_cpd, lBand_mu, '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.2, 'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Late (220-355 ms)');
yline(eNum, '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.2, 'DisplayName', 'Number (Early ref)');
yline(lNum, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, 'DisplayName', 'Number (Late ref)');
xline(optBandSF_cpd, 'r:', 'LineWidth', 1.8, 'DisplayName', sprintf('Ideal Observer (%.2f cpd)', optBandSF_cpd));
set(gca, 'XScale', 'log'); xlim([sfFreqs_cpd(1)*0.8, sfFreqs_cpd(end)*1.25]);
xticks(sfFreqs_cpd); xticklabels(arrayfun(@(f) sprintf('%.2f', f), sfFreqs_cpd, 'UniformOutput', false));
xlabel('Spatial Frequency (cpd)'); ylabel('Mean Incremental \DeltaR^2');
title('Spatial-Frequency Tuning: Early vs. Late Window');
legend('Location', 'northwest'); grid on;
save_fig_png(fS5, fullfile(suppDir, 'supp_5_sf_tuning_early_vs_late.png'), 300);
close(fS5);

fprintf('Saved 4 main panels to:         %s\n', mainDir);
fprintf('Saved 5 supplementary panels to: %s\n', suppDir);

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
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
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
