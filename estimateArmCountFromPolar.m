function armCount = estimateArmCountFromPolar(polarImg)
% Improved estimate of the number of spiral arms from the polar image
% - Higher MinPeakHeight to reduce noise
% - Ignore outlier rings
% - Easy to tune parameters

numRadii = 12; % Number of radial samples (rings)
minPeakHeight = 0.32; % Increased from 0.25 for robustness
maxArms = 8; % Ignore rings with too many peaks
minArms = 2; % Ignore rings with only 1 peak (often noise)
[rows, cols] = size(polarImg);
radialIndices = round(linspace(round(rows*0.3), round(rows*0.85), numRadii));
peakCounts = zeros(numRadii, 1);

for i = 1:numRadii
    ring = polarImg(radialIndices(i), :);
    ring = ring - min(ring); % Remove background
    ring = ring / max(ring + eps); % Normalize
    % Smooth to reduce noise
    ring = smoothdata(ring, 'gaussian', 7);
    % Find peaks above a threshold
    [pks, locs] = findpeaks(ring, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', round(cols/10));
    peakCounts(i) = numel(pks);
end

% Remove outlier rings (too few or too many peaks)
validCounts = peakCounts(peakCounts >= minArms & peakCounts <= maxArms);

% Use the mode of the counts (most common value) as the arm count
if ~isempty(validCounts)
    armCount = mode(validCounts);
else
    armCount = 0; % Cannot determine
end
end
