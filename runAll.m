%% RUNALL
% Master pipeline script:
%   1) Auto-detects any newly completed MEG .7z archives (>= 1 GB, ignoring .part)
%      and preprocesses only new subjects into 32x32 gradiometer RDMs.
%   2) Derives the ideal-observer optimal SF (if not already cached) from
%      CSF-weighted Fourier power differences across all numerosity pairs (>0.5 cpd).
%   3) Builds the 5-band Contrast Energy model RDMs with Band 3 aligned to optimalSF.
%   4) Runs the group MEG replication, variance partitioning, and SF analysis.

clear; close all; clc;

cd(fileparts(mfilename('fullpath')));

megDir    = fullfile(pwd, 'MEG_Data');
megRdmDir = fullfile(pwd, 'meg_rdms');
pyScript  = fullfile(pwd, 'preprocess_meg_to_rdms.py');

%% 1. Auto-detect complete .7z archives and preprocess any new subjects
archives         = dir(fullfile(megDir, 'sub-*.7z'));
isComplete       = [archives.bytes] >= 1e9; % strictly complete archives (>= 1 GB)
completeArchives = archives(isComplete);
partFiles        = dir(fullfile(megDir, 'sub-*.7z.part'));

fprintf('=== Step 1: Checking MEG_Data folder ===\n');
fprintf('  Complete .7z archives found: %d | In-progress .part downloads: %d\n', ...
    numel(completeArchives), numel(partFiles));

newSubjects = {};
for k = 1:numel(completeArchives)
    [~, subjName] = fileparts(completeArchives(k).name);
    rdmFile = fullfile(megRdmDir, [subjName '_rdm_grad.mat']);
    if ~exist(rdmFile, 'file')
        newSubjects{end+1} = subjName; %#ok<SAGROW>
    end
end

if ~isempty(newSubjects)
    fprintf('  New subject(s) to preprocess: %s\n', strjoin(newSubjects, ', '));
    cmd = sprintf('python3 "%s" --meg-dir "%s" --out-dir "%s" --subjects %s', ...
        pyScript, megDir, megRdmDir, strjoin(newSubjects, ' '));
    fprintf('  Running MEG preprocessing via Python/MNE...\n\n');
    status = system(cmd, '-echo');
    assert(status == 0, 'preprocess_meg_to_rdms.py exited with status %d.', status);
else
    fprintf('  All %d complete subject(s) are already preprocessed in meg_rdms/.\n', numel(completeArchives));
end

%% 2. Ideal-observer optimal SF derivation (auto-recomputes if CSF or stim params changed)
fprintf('\n=== Step 2: Checking ideal-observer optimal SF ===\n');
paramsNow = frf_default_params();
optFile   = fullfile(pwd, 'results', 'optimal_sf.mat');
needRecompute = true;
if exist(optFile, 'file')
    L = load(optFile, 'opt');
    if isfield(L, 'opt') && isfield(L.opt, 'csfParams') && isfield(L.opt, 'stimParams') ...
            && isequaln(L.opt.csfParams, paramsNow.csf) && isequaln(L.opt.stimParams, paramsNow.stim)
        needRecompute = false;
        fprintf('  Using cached optimal SF = %.3f cpd (matches current CSF & stim params).\n', ...
            L.opt.optimalSF_cpd);
    else
        fprintf('  CSF or stimulus parameters changed; recomputing optimal SF...\n');
    end
end
if needRecompute
    derive_optimal_sf(paramsNow);
end

%% 3. Build Contrast Energy model RDMs (center band = optimal SF)
fprintf('\n=== Step 3: Running Contrast Energy model simulation ===\n');
run_frf_rdm_simulation;

%% 4. Run group MEG replication, variance partitioning, and SF analysis
fprintf('\n=== Step 4: Running MEG group analysis across all processed subjects ===\n');
run_meg_vp_test;

%% 5. Regenerate ideal-observer walkthrough figures
fprintf('\n=== Step 5: Regenerating ideal-observer walkthrough figures ===\n');
plot_ideal_observer_walkthrough;

fprintf('\n=== runAll complete! ===\n');
