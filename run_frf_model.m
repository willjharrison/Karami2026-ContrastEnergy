function out = run_frf_model(im, params, precomp)
% RUN_FRF_MODEL Single-stage local contrast-energy model built directly on
% Will Harrison's createBandPassFilters_customSizes.m and polarDistFun.m.
%
% Tiles spatial-frequency space into evenly spaced (1-octave) raised-log-cosine
% bands, rectifies each filtered image to compute local contrast energy, and
% spatially pools each band.
%
% Inputs:
%   im      - Square 2D double image (imSize x imSize)
%   params  - Parameter struct from frf_default_params()
%   precomp - (Optional) Precomputed filter bank from build_frf_filters(params)
%
% Output (struct):
%   out.cePooled    - [1 x nBands] spatially pooled contrast energy per SF band
%   out.sfFreqs_cpd - [1 x nBands] center frequencies in cycles/degree
%   out.ceMaps      - [imSize x imSize x nBands] (only if params.returnMaps == true)

if nargin < 3 || isempty(precomp)
    precomp = build_frf_filters(params);
end

returnMaps = isfield(params, 'returnMaps') && params.returnMaps;

imIn = double(im);
if params.sf.subtractMean
    imIn = imIn - mean(imIn(:));
end
imFFT = fft2(imIn);

useBands = params.sf.useBands;
nBands   = numel(useBands);
imSize   = size(imIn, 1);

cePooled = zeros(1, nBands);
if returnMaps
    ceMaps = zeros(imSize, imSize, nBands);
end

for i = 1:nBands
    fShift = precomp.sfFiltsShift(:, :, useBands(i));
    resp   = real(ifft2(imFFT .* fShift));
    enMap  = apply_rectify(resp, params.sf.rectifyMode);
    
    cePooled(i) = pool_map(enMap, params.sf.poolMode);
    if returnMaps
        ceMaps(:, :, i) = enMap;
    end
end

out.cePooled    = cePooled;
out.sfFreqs_cpd = precomp.sfFreqs_cpd(useBands);
if returnMaps
    out.ceMaps = ceMaps;
end
end

function y = apply_rectify(x, mode)
switch lower(mode)
    case 'abs',      y = abs(x);
    case 'square',   y = x .^ 2;
    case 'halfwave', y = max(x, 0);
    otherwise,       error('Unknown rectifyMode: %s', mode);
end
end

function v = pool_map(m, mode)
switch lower(mode)
    case 'mean',     v = mean(m(:));
    case 'sum',      v = sum(m(:));
    case 'rms',      v = sqrt(mean(m(:).^2));
    case 'log_mean', v = log(mean(m(:)) + eps);
    otherwise,       error('Unknown poolMode: %s', mode);
end
end
