function params = frf_default_params()
% FRF_DEFAULT_PARAMS Default parameter configuration for the single-stage
% Contrast Energy (CE) model applied to Karami et al. (2026, Nature Comms) displays.
%
% Uses an ideal-observer, CSF-weighted Fourier power simulation (derive_optimal_sf.m)
% to set the center spatial frequency (Band 3 of 5) a priori without touching MEG data.

%% Paths (relative to repository root)
repoDir               = fileparts(mfilename('fullpath'));
params.filterDir      = fullfile(repoDir, 'filters');
params.outputDir      = fullfile(repoDir, 'results');

%% Karami et al. (2026) stimulus grid (32 conditions: 4 N x 4 AIA x 2 TFA)
% Condition order matches meg_mne_preprocess.ipynb: N1S1TFA1, N1S1TFA2, ..., N4S4TFA2
params.stim.numerosities    = [6, 10, 17, 29];
params.stim.itemAreas_vd2   = [0.04, 0.07, 0.12, 0.20]; % visual deg^2
params.stim.fieldDiams_deg  = [9.0, 13.5];              % virtual circle diameter (deg)

% Canvas geometry (square power-of-2 for fast FFT and createBandPassFilters_customSizes)
params.stim.imSize_px       = 512;   % pixels (square)
params.stim.imSize_deg      = 16.0;  % deg of visual angle (> 13.5 deg to pad aperture edges)
params.stim.ppd             = params.stim.imSize_px / params.stim.imSize_deg; % 32 px/deg

% Dot rendering
params.stim.bgLum           = 0.5;   % mid-gray background
params.stim.dotLum          = 0.0;   % black dots
params.stim.dotAreaJitter   = 0.50;  % +/-50% within-array dot-area jitter (mean exact = itemArea_vd2)
params.stim.minEdgeSep_deg  = 0.10;  % minimum edge-to-edge gap between dots (deg)
params.stim.includeFixCross = false; % red cross is constant across conditions
params.stim.nExemplars      = 100;   % random dot arrays averaged per condition
params.stim.rngSeed         = 42;    % reproducible synthetic stimuli

%% Uncontroversial Human Photopic CSF (createCSF.m: truncated log-parabola)
% Standard photopic foveal/parafoveal sensitivity (peak ~3 cpd, ~2-oct full bandwidth)
params.csf.fmax             = 3.0;   % peak spatial frequency (cpd)
params.csf.vmax             = 1.0;   % normalized peak gain = 1
params.csf.beta             = 1.0;   % 1-octave half-bandwidth parameter
params.csf.trunc            = 0.85;  % low-SF truncation (floor = 0.15)
params.csf.minSF_cpd        = 0.5;   % lower SF bound for ideal-observer peak search
params.csf.nDisplaysPerNum  = 200;   % displays per numerosity in derive_optimal_sf.m

%% Spatial-Frequency Filter Bank (5 bands, Band 3 aligned to optimal SF)
params.sf.useOptimalSF      = true;  % align Band 3 to ideal-observer optimalSF_cpd
params.sf.nFreqs            = 5;     % 5 isotropic bands
params.sf.stepOctaves       = 1.0;   % 1-octave spacing around Band 3
params.sf.lowPassFreq_cpd   = 0.5;   % fallback if useOptimalSF == false
params.sf.highPassFreq_cpd  = 8.0;   % fallback if useOptimalSF == false
params.sf.subtractMean      = true;  % remove image DC before FFT
params.sf.zeroDC_lp         = true;  % zero DC bin of lpFilter
params.sf.useBands          = 1:params.sf.nFreqs;

% Rectification & spatial pooling (simplest contrast-energy computation)
params.sf.rectifyMode       = 'abs';
params.sf.poolMode          = 'mean';

%% RDM construction
params.rdm.distMetric       = 'log_euclidean';
params.rdm.zscoreFeatures   = false;

%% Diagnostics & Plotting
params.plotSanityChecks     = true;
params.saveResults          = true;
end
