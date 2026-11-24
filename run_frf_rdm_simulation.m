%% RUN_FRF_RDM_SIMULATION
% Builds 32x32 Contrast Energy model RDMs for the Karami et al. (2026)
% stimulus grid.
%   - Primary model (rdm.ContrastEnergy): single spatial-frequency band
%     aligned to the CSF-limited ideal-observer peak from derive_optimal_sf.m.
%   - Supplementary models (rdm.SF_Bands, rdm.ContrastEnergy_MultiBand):
%     full multi-band filter bank spanning surrounding spatial frequencies.
%
% Outputs (saved to FRF_model/results/):
%   - fig1_filter_tiling.png, fig2_stimuli_and_ce_maps.png, fig3_model_rdms.png
%   - rdm_correlations.csv, summary.txt, and frf_rdms.mat

clear; close all; clc;
tStart = tic;

%% 1. Load parameters and align 5-band SF bank to ideal-observer optimal SF
params = frf_default_params();
addpath(params.filterDir);
if params.saveResults && ~exist(params.outputDir, 'dir')
    mkdir(params.outputDir);
end

nyquist_cpd = (params.stim.imSize_px / 2) / params.stim.imSize_deg;
if params.sf.useOptimalSF
    optFile = fullfile(params.outputDir, 'optimal_sf.mat');
    needRecompute = true;
    if exist(optFile, 'file')
        L = load(optFile, 'opt');
        if isfield(L, 'opt') && isfield(L.opt, 'csfParams') && isfield(L.opt, 'stimParams') ...
                && isequaln(L.opt.csfParams, params.csf) && isequaln(L.opt.stimParams, params.stim)
            opt = L.opt;
            needRecompute = false;
        end
    end
    if needRecompute
        opt = derive_optimal_sf(params);
    end
    targetSF = opt.optimalSF_cpd;
    
    % Center middle band of nFreqs on targetSF with params.sf.stepOctaves spacing
    % (clamping octave step slightly only if highest band would exceed Nyquist)
    halfSpan = (params.sf.nFreqs - 1) / 2;
    maxOct   = floor(log2(0.95 * nyquist_cpd / targetSF) / halfSpan * 100) / 100;
    stepOct  = min(params.sf.stepOctaves, maxOct);
    params.sf.lowPassFreq_cpd  = targetSF / (2 ^ (halfSpan * stepOct));
    params.sf.highPassFreq_cpd = targetSF * (2 ^ (halfSpan * stepOct));
    fprintf('Aligned %d-band SF bank: Band %d = %.3f cpd (step = %.2f octaves)\n', ...
        params.sf.nFreqs, ceil(params.sf.nFreqs / 2), targetSF, stepOct);
end

precomp = build_frf_filters(params);
sfFreqs_cpd = precomp.sfFreqs_cpd(params.sf.useBands);
octaveStep  = log2(sfFreqs_cpd(2) / sfFreqs_cpd(1));

fprintf('=== Contrast Energy Model Simulation (Karami et al. 32-condition grid) ===\n');
fprintf('SF Bank: %d bands (%.3f - %.2f cpd, step = %.2f octaves), rectify=%s, pool=%s, dist=%s\n', ...
    numel(sfFreqs_cpd), sfFreqs_cpd(1), sfFreqs_cpd(end), octaveStep, ...
    params.sf.rectifyMode, params.sf.poolMode, params.rdm.distMetric);
fprintf('Band centers (cpd): %s\n', mat2str(round(sfFreqs_cpd, 3)));

%% 2. Build 32-condition table in exact Karami et al. order (N -> AIA -> TFA)
nums  = params.stim.numerosities;
areas = params.stim.itemAreas_vd2;
diams = params.stim.fieldDiams_deg;

nCond = numel(nums) * numel(areas) * numel(diams);
condTable = table('Size', [nCond, 7], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Label', 'Number', 'ItemArea', 'FieldDiam', 'FieldArea', 'SurfaceArea', 'Density'});

cIdx = 0;
for iN = 1:numel(nums)
    for iS = 1:numel(areas)
        for iF = 1:numel(diams)
            cIdx = cIdx + 1;
            N   = nums(iN);
            AIA = areas(iS);
            FD  = diams(iF);
            TFA = pi * (FD / 2)^2;
            condTable.Label(cIdx)       = sprintf('N%dS%dTFA%d', iN, iS, iF);
            condTable.Number(cIdx)      = N;
            condTable.ItemArea(cIdx)    = AIA;
            condTable.FieldDiam(cIdx)   = FD;
            condTable.FieldArea(cIdx)   = TFA;
            condTable.SurfaceArea(cIdx) = N * AIA;
            condTable.Density(cIdx)     = N / TFA;
        end
    end
end

%% 3. Visual sanity checks on filter bank & example stimulus
rng(params.stim.rngSeed);
if params.plotSanityChecks
    plot_sanity_checks(params, precomp, condTable);
end

%% 4. Compute contrast energy across 32 conditions x nExemplars
nBands = numel(params.sf.useBands);
ceCondMean  = zeros(nCond, nBands);
lumCondMean = zeros(nCond, 1);

fprintf('\nRunning Contrast Energy bank across %d conditions x %d exemplars (%d images)...\n', ...
    nCond, params.stim.nExemplars, nCond * params.stim.nExemplars);

for c = 1:nCond
    ceAccum  = zeros(params.stim.nExemplars, nBands);
    lumAccum = zeros(params.stim.nExemplars, 1);
    for e = 1:params.stim.nExemplars
        im = generate_karami_stimulus( ...
            condTable.Number(c), condTable.ItemArea(c), condTable.FieldDiam(c), params.stim);
        res = run_frf_model(im, params, precomp);
        ceAccum(e, :) = res.cePooled;
        lumAccum(e)   = mean(im(:));
    end
    ceCondMean(c, :) = mean(ceAccum, 1);
    lumCondMean(c)   = mean(lumAccum);
    
    if mod(c, 8) == 0 || c == nCond
        fprintf('  Condition %2d/%d (%s) done [%.1f s elapsed]\n', ...
            c, nCond, condTable.Label(c), toc(tStart));
    end
end

%% 5. Construct 32x32 Model RDMs
rdm.Number         = abs(log(condTable.Number)      - log(condTable.Number)');
rdm.ItemArea       = abs(log(condTable.ItemArea)    - log(condTable.ItemArea)');
rdm.FieldArea      = abs(log(condTable.FieldArea)   - log(condTable.FieldArea)');
rdm.SurfaceArea    = abs(log(condTable.SurfaceArea) - log(condTable.SurfaceArea)');
rdm.Density        = abs(log(condTable.Density)     - log(condTable.Density)');

contourProxy       = condTable.Number .* sqrt(condTable.ItemArea);
rdm.Contour        = abs(log(contourProxy) - log(contourProxy)');

rdm.SF_Bands = cell(1, nBands);
sfBandNames  = cell(1, nBands);
for b = 1:nBands
    v = log(max(ceCondMean(:, b), eps));
    rdm.SF_Bands{b} = abs(v - v');
    sfBandNames{b}  = sprintf('SF_%.2fcpd', sfFreqs_cpd(b));
end

% Primary model: single spatial-frequency band aligned to ideal-observer peak
if exist('targetSF', 'var')
    [~, optBandIdx] = min(abs(sfFreqs_cpd - targetSF));
else
    optBandIdx = ceil(nBands / 2);
end
optBandSF_cpd = sfFreqs_cpd(optBandIdx);

rdm.ContrastEnergy           = rdm.SF_Bands{optBandIdx};
rdm.ContrastEnergy_MultiBand = compute_feature_rdm(ceCondMean, params.rdm);
rdm.MeanLuminance            = abs(lumCondMean - lumCondMean');

%% 6. Collinearity & Partial-r Diagnostics vs. Number RDM
baseNames = {'Number', 'ItemArea', 'FieldArea', 'SurfaceArea', 'Density', 'Contour', 'ContrastEnergy', 'ContrastEnergy_MultiBand', 'MeanLuminance'};
utMask = triu(true(nCond), 1);
nPairs = sum(utMask(:));

RDM_base = zeros(nPairs, numel(baseNames));
for m = 1:numel(baseNames)
    M = rdm.(baseNames{m});
    RDM_base(:, m) = M(utMask);
end

corrMat = corr(RDM_base, 'Type', 'Pearson');
corrTbl = array2table(round(corrMat, 3), 'RowNames', baseNames, 'VariableNames', baseNames);

fprintf('\n=== Core RDM Pearson Correlation Matrix (496 condition pairs) ===\n');
disp(corrTbl);

X_nonNum = RDM_base(:, 2:5);
y_num    = RDM_base(:, 1);
resNum   = y_num - [ones(nPairs,1), X_nonNum] * ([ones(nPairs,1), X_nonNum] \ y_num);

summaryLines = {};
summaryLines{end+1} = sprintf('Primary Model: Single ideal-observer SF band = Band %d (%.3f cpd)', optBandIdx, optBandSF_cpd);
summaryLines{end+1} = sprintf('Supplementary SF Bank: %d bands (%.3f - %.2f cpd, %.2f-octave steps), rectify=%s, pool=%s, dist=%s', ...
    nBands, sfFreqs_cpd(1), sfFreqs_cpd(end), octaveStep, ...
    params.sf.rectifyMode, params.sf.poolMode, params.rdm.distMetric);
summaryLines{end+1} = '';
summaryLines{end+1} = '=== Model Diagnostics vs. Number RDM ===';

fprintf('=== Model Diagnostics vs. Number RDM ===\n');
for idx = [6, 7, 8, 9]
    r0 = corr(y_num, RDM_base(:, idx));
    rc = RDM_base(:, idx) - [ones(nPairs,1), X_nonNum] * ([ones(nPairs,1), X_nonNum] \ RDM_base(:, idx));
    rp = corr(resNum, rc);
    line = sprintf('  %-15s | zero-order r(Num) = %+.3f | partial r|non-num = %+.3f', baseNames{idx}, r0, rp);
    fprintf('%s\n', line);
    summaryLines{end+1} = line;
end

fprintf('\n=== Individual SF Bands vs. Number RDM ===\n');
summaryLines{end+1} = '';
summaryLines{end+1} = '=== Individual SF Bands vs. Number RDM ===';
for b = 1:nBands
    vb = rdm.SF_Bands{b}(utMask);
    r0 = corr(y_num, vb);
    rc = vb - [ones(nPairs,1), X_nonNum] * ([ones(nPairs,1), X_nonNum] \ vb);
    rp = corr(resNum, rc);
    line = sprintf('  Band %d (%6.2f cpd) | zero-order r(Num) = %+.3f | partial r|non-num = %+.3f', ...
        b, sfFreqs_cpd(b), r0, rp);
    fprintf('%s\n', line);
    summaryLines{end+1} = line;
end

fprintf('\nTotal runtime: %.1f s\n', toc(tStart));

%% 7. Plot RDM comparison (saved silently to results/panels_supplement/)
if params.plotSanityChecks
    suppDir = fullfile(params.outputDir, 'panels_supplement');
    if ~exist(suppDir, 'dir'), mkdir(suppDir); end
    fig3 = figure('Visible', 'off', 'Name', '32x32 Model RDMs', 'Color', 'w', 'Position', [100 100 1200 650]);
    showModels = {'Number', 'FieldArea', 'SurfaceArea', 'ItemArea', 'Contour', 'ContrastEnergy'};
    showTitles = {'Number', 'FieldArea', 'SurfaceArea', 'ItemArea', 'Contour', sprintf('Contrast Energy (%.2f cpd)', optBandSF_cpd)};
    for k = 1:numel(showModels)
        subplot(2, 3, k);
        imagesc(rdm.(showModels{k}));
        axis square; colorbar;
        title(showTitles{k});
        xlabel('Condition'); ylabel('Condition');
    end
    if params.saveResults
        save_fig_png(fig3, fullfile(suppDir, 'supp_model_rdms.png'));
    end
    close(fig3);
end

%% 8. Save results
if params.saveResults
    saveFile = fullfile(params.outputDir, 'frf_rdms.mat');
    save(saveFile, 'rdm', 'condTable', 'ceCondMean', 'lumCondMean', 'sfFreqs_cpd', 'sfBandNames', 'optBandIdx', 'optBandSF_cpd', 'params', 'corrTbl');
    writetable(corrTbl, fullfile(params.outputDir, 'rdm_correlations.csv'), 'WriteRowNames', true);
    fid = fopen(fullfile(params.outputDir, 'summary.txt'), 'w');
    if fid ~= -1
        fprintf(fid, '%s\n', summaryLines{:});
        fclose(fid);
    end
    fprintf('Saved Contrast Energy RDMs, CSV, summary.txt, and PNG figures to: %s\n', params.outputDir);
end

%% Helper functions
function D = compute_feature_rdm(X, rdmParams)
if rdmParams.zscoreFeatures
    X = zscore(X, 0, 1);
end
switch lower(rdmParams.distMetric)
    case 'correlation'
        D = 1 - corr(X');
        D(1:size(D,1)+1:end) = 0;
    case 'euclidean'
        D = squareform(pdist(X, 'euclidean'));
    case 'log_euclidean'
        D = squareform(pdist(log(max(X, eps)), 'euclidean'));
    otherwise
        error('Unknown rdm.distMetric: %s', rdmParams.distMetric);
end
end

function plot_sanity_checks(params, precomp, condTable)
imSize = params.stim.imSize_px;
midRow = imSize / 2 + 1;
xCols  = (imSize / 2 + 1):imSize;
fAxis  = precomp.radDist_cpd(midRow, xCols);
nBands = size(precomp.sfFilts, 3);
cols   = parula(nBands + 1);

suppDir = fullfile(params.outputDir, 'panels_supplement');
if ~exist(suppDir, 'dir'), mkdir(suppDir); end

fig1 = figure('Visible', 'off', 'Name', 'SF Filter Bank Tiling', 'Color', 'w', 'Position', [100 420 750 420]);
hold on;
for k = 1:nBands
    plot(fAxis, precomp.sfFilts(midRow, xCols, k), 'Color', cols(k, :), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Band %d: %.2f cpd', k, precomp.sfFreqs_cpd(k)));
end
plot(fAxis, sum(precomp.sfFilts(midRow, xCols, :), 3), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Partition Sum');
set(gca, 'XScale', 'log'); xlim([fAxis(2), max(fAxis)]); ylim([-0.05 1.15]);
xlabel('Spatial Frequency (cpd)'); ylabel('Filter Gain');
title(sprintf('Raised-Log-Cosine Filter Bank (n=%d, Center = %.2f cpd)', nBands, precomp.sfFreqs_cpd(ceil(nBands/2))));
legend('Location', 'eastoutside'); grid on;
if params.saveResults
    save_fig_png(fig1, fullfile(suppDir, 'supp_filter_tiling.png'));
end
close(fig1);

ci = 32; % N4S4TFA2 (29 large dots, 13.5 deg)
im = generate_karami_stimulus(condTable.Number(ci), condTable.ItemArea(ci), condTable.FieldDiam(ci), params.stim);
pMaps = params; pMaps.returnMaps = true;
exOut = run_frf_model(im, pMaps, precomp);

fig2 = figure('Visible', 'off', 'Name', 'Example Stimulus & Contrast Energy Maps', 'Color', 'w', 'Position', [100 80 1100 500]);
subplot(2, 3, 1); imshow(im, [0 1]); title('Input (N=29, A=0.20, D=13.5)');
for k = 1:min(5, nBands)
    subplot(2, 3, k + 1);
    imagesc(exOut.ceMaps(:, :, k)); axis image off;
    title(sprintf('Band %d (%.2f cpd)', k, exOut.sfFreqs_cpd(k)));
end
if params.saveResults
    save_fig_png(fig2, fullfile(suppDir, 'supp_stimuli_ce_maps.png'));
end
close(fig2);
end

function save_fig_png(figHandle, filepath)
if exist('exportgraphics', 'file')
    exportgraphics(figHandle, filepath, 'Resolution', 150);
else
    saveas(figHandle, filepath);
end
end
