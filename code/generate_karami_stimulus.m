function [im, meta] = generate_karami_stimulus(nDots, itemArea_vd2, fieldDiam_deg, stimParams)
% GENERATE_KARAMI_STIMULUS Generates a synthetic dot array matching
% Karami et al. (2026, Nature Comms) stimulus specifications, including
% within-array dot-size heterogeneity constrained to exact condition-mean
% Average Item Area (and therefore exact Total Surface Area).
%
% Inputs:
%   nDots         - integer numerosity (6, 10, 17, or 29)
%   itemArea_vd2  - average individual dot area in visual deg^2
%   fieldDiam_deg - virtual circular aperture diameter in visual deg
%   stimParams    - struct from frf_default_params().stim
%
% Outputs:
%   im   - imSize_px x imSize_px double image (bgLum background, dotLum dots)
%   meta - struct with actual placed dot centers (px) and radii (px)

imSize = stimParams.imSize_px;
ppd    = stimParams.ppd;

% Sample heterogeneous dot areas around itemArea_vd2 (preserving exact mean)
if isfield(stimParams, 'dotAreaJitter') && stimParams.dotAreaJitter > 0
    j = stimParams.dotAreaJitter;
    rawWeights = (1 - j) + (2 * j) * rand(nDots, 1);
    dotAreas_vd2 = itemArea_vd2 * (rawWeights / mean(rawWeights));
else
    dotAreas_vd2 = ones(nDots, 1) * itemArea_vd2;
end

% Convert degrees to pixels; sort largest dots first for reliable packing
dotRadii_deg   = sqrt(dotAreas_vd2 / pi);
dotRadii_px    = sort(dotRadii_deg * ppd, 'descend');
fieldRadius_px = (fieldDiam_deg / 2) * ppd;
minEdgeSep_px  = stimParams.minEdgeSep_deg * ppd;

% Place non-overlapping dots inside the virtual circular aperture
centers     = zeros(nDots, 2);
placedRadii = zeros(nDots, 1);
placed      = 0;
maxAttempts = 15000;
attempts    = 0;

while placed < nDots && attempts < maxAttempts
    attempts = attempts + 1;
    rDot = dotRadii_px(placed + 1);
    maxRadiusPlace = max(0, fieldRadius_px - rDot);
    
    % Uniform sampling over disk area
    r = maxRadiusPlace * sqrt(rand());
    theta = 2 * pi * rand();
    cand = [r * cos(theta), r * sin(theta)];
    
    if placed == 0
        placed = 1;
        centers(1, :)  = cand;
        placedRadii(1) = rDot;
    else
        dists   = sqrt(sum((centers(1:placed, :) - cand).^2, 2));
        minSeps = placedRadii(1:placed) + rDot + minEdgeSep_px;
        if all(dists >= minSeps)
            placed = placed + 1;
            centers(placed, :)  = cand;
            placedRadii(placed) = rDot;
        end
    end
end

if placed < nDots
    warning('generate_karami_stimulus:PlacementIncomplete', ...
        'Placed %d/%d dots (N=%d, area=%.2f, field=%.1f). Reduce minEdgeSep_deg if needed.', ...
        placed, nDots, nDots, itemArea_vd2, fieldDiam_deg);
    centers     = centers(1:placed, :);
    placedRadii = placedRadii(1:placed);
end

% Render image
im = ones(imSize, imSize) * stimParams.bgLum;
halfSz = imSize / 2;
[X, Y] = meshgrid((1:imSize) - halfSz - 0.5, (1:imSize) - halfSz - 0.5);

for d = 1:size(centers, 1)
    mask = ((X - centers(d, 1)).^2 + (Y - centers(d, 2)).^2) <= placedRadii(d)^2;
    im(mask) = stimParams.dotLum;
end

% Optional thin fixation cross (luminance proxy; in the actual MEG study it was red on gray)
if isfield(stimParams, 'includeFixCross') && stimParams.includeFixCross
    armLen_px   = round(0.5 * ppd);
    armWidth_px = max(1, round(0.05 * ppd));
    crossMask   = (abs(X) <= armLen_px & abs(Y) <= armWidth_px) | ...
                  (abs(Y) <= armLen_px & abs(X) <= armWidth_px);
    im(crossMask) = 0.75; % placeholder contrast for luminance-channel test
end

meta.nPlaced        = size(centers, 1);
meta.centers_px     = centers;
meta.dotRadii_px    = placedRadii;
meta.dotRadius_px   = mean(placedRadii);
meta.fieldRadius_px = fieldRadius_px;
end
