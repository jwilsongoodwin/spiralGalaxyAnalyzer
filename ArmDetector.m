clear all
close all
function analyzeGalaxySpiralManualCut()
    % Select a galaxy image
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.tif;*.bmp', 'Images (*.jpg, *.png, *.tif, *.bmp)'}, 'Select a Galaxy Image');
    if isequal(filename, 0)
        disp('User cancelled.');
        return;
    end
    fullPath = fullfile(pathname, filename);
    img = imread(fullPath);

    % Convert to grayscale if needed
    if size(img, 3) > 1
        grayImg = rgb2gray(img);
    else
        grayImg = img;
    end
    grayImg = im2double(grayImg);

    % Resize if too large
    if max(size(grayImg)) > 1200
        grayImg = imresize(grayImg, 1200/max(size(grayImg)));
    end

    % Logarithmic stretch
    grayImg = log1p(grayImg);
    grayImg = (grayImg - min(grayImg(:))) / (max(grayImg(:)) - min(grayImg(:)));

    % Show image and let user pick center manually
    figure('Name', 'Select Galaxy Center and Mask Core');
    imshow(grayImg, []);
    title('Click on the center of the galaxy, then adjust circle to mask core');
    [centerX, centerY] = ginput(1);

    % Draw initial circle mask with default radius
    h = drawcircle('Center', [centerX, centerY], 'Radius', 30, 'Color', 'r');
    disp('Adjust the circle as needed, then double-click inside to confirm.');
    wait(h);
    mask = createMask(h);
    maskedImg = grayImg;
    maskedImg(mask) = 0; % Mask out core

    % Preview the masked image
    figure('Name', 'Preview Masked Image');
    imshow(maskedImg, []);
    title('Masked Image (Core Removed)');
    pause(2);
    close;

    % Make image square around center
    [rows, cols] = size(maskedImg);
    maxDim = max(rows, cols);
    squareImg = zeros(maxDim, maxDim);
    rowOffset = floor((maxDim - rows) / 2);
    colOffset = floor((maxDim - cols) / 2);
    squareImg((1:rows) + rowOffset, (1:cols) + colOffset) = maskedImg;
    centerX = centerX + colOffset;
    centerY = centerY + rowOffset;

    % Polar transform AFTER masking
    maxRadius = min([centerX, centerY, maxDim - centerX, maxDim - centerY]);
    [polarImg, ~, ~] = cartesian2Polar(squareImg, centerX, centerY, maxRadius);

    % Normalize polar image before FFT
    polarImg = (polarImg - min(polarImg(:))) / (max(polarImg(:)) - min(polarImg(:)));

    % Preview polar image before FFT
    figure('Name', 'Polar Transformed Image');
    imshow(polarImg, []);
    title('Polar Image (Before FFT)');
    pause(2);
    close;

    % 2D FFT of polar image
    fftPolar = fft2(polarImg);
    fftShift = fftshift(fftPolar);

    % High-pass filter in Fourier space
    [rows, cols] = size(fftShift);
    [X, Y] = meshgrid(1:cols, 1:rows);
    centerRow = ceil(rows/2);
    centerCol = ceil(cols/2);
    hpRadius = 3; % Reduced high-pass radius
    hpMask = sqrt((X - centerCol).^2 + (Y - centerRow).^2) > hpRadius;
    fftShift = fftShift .* hpMask;

    % Analyze spiral
    [pitchAngle, barStrength, dominantMode] = analyzeImprovedSpiralFFT(fftShift, polarImg);

    % Display results
    fprintf('Dominant Spiral Mode (Number of Arms): %d\n', dominantMode);
    fprintf('Estimated Pitch Angle: %.2f degrees\n', pitchAngle);
    fprintf('Bar Strength: %.2f\n', barStrength);
end

function [polarImg, rVals, thetaVals] = cartesian2Polar(img, centerX, centerY, maxRadius)
    numRadialBins = 200;
    numAngularBins = 360;
    rVals = linspace(1, maxRadius, numRadialBins);
    thetaVals = linspace(0, 2*pi, numAngularBins);

    [THETA, R] = meshgrid(thetaVals, rVals);
    X = centerX + R .* cos(THETA);
    Y = centerY + R .* sin(THETA);

    polarImg = interp2(img, X, Y, 'linear', 0);
end

function [pitchAngle, barStrength, dominantMode] = analyzeImprovedSpiralFFT(fftData, polarImg)
    [rows, cols] = size(fftData);
    centerRow = ceil(rows/2);
    centerCol = ceil(cols/2);

    % --- Dominant Mode Detection (Improved) ---
    angularProfile = abs(fftData(centerRow, :));
    angularProfile = movmean(angularProfile, 5);

    mValues = 1:10;
    powerByMode = zeros(size(mValues));
    for i = 1:length(mValues)
        offset = mValues(i);
        leftIdx = centerCol - offset;
        rightIdx = centerCol + offset;
        if leftIdx > 0 && rightIdx <= cols
            power = angularProfile(leftIdx) + angularProfile(rightIdx);
            powerByMode(i) = power / mValues(i);
        end
    end

    [~, maxIdx] = max(powerByMode);
    dominantMode = mValues(maxIdx);

    % --- Pitch Angle Calculation ---
    phase = angle(fftData);
    modeCol = centerCol + dominantMode;
    if modeCol > cols
        modeCol = modeCol - cols;
    end
    radialPhaseProfile = unwrap(phase(:, modeCol));
    midPoint = floor(length(radialPhaseProfile) / 2);
    radialRange = max(1, midPoint-20):min(length(radialPhaseProfile), midPoint+20);
    x = radialRange';
    y = radialPhaseProfile(radialRange);
    p = polyfit(x, y, 1);
    slope = p(1);

    pitchAngle1 = atand(abs(dominantMode / (slope * rows/4)));

    % --- Spatial Domain (Log-polar Hough Transform) ---
    [polarRows, polarCols] = size(polarImg);
    edges = edge(polarImg, 'canny');
    [H, theta, rho] = hough(edges);
    peaks = houghpeaks(H, 5, 'Threshold', 0.5*max(H(:)));
    lines = houghlines(edges, theta, rho, peaks, 'FillGap', 20, 'MinLength', 30);

    if ~isempty(lines)
        lineAngles = zeros(length(lines), 1);
        for i = 1:length(lines)
            dx = lines(i).point2(1) - lines(i).point1(1);
            dy = lines(i).point2(2) - lines(i).point1(2);
            lineAngles(i) = atand(dy/dx);
        end
        validAngles = lineAngles(lineAngles > 5 & lineAngles < 85);
        if ~isempty(validAngles)
            pitchAngle2 = mean(validAngles);
        else
            pitchAngle2 = pitchAngle1;
        end
    else
        pitchAngle2 = pitchAngle1;
    end

    pitchAngle = 0.7 * pitchAngle1 + 0.3 * pitchAngle2;
    pitchAngle = min(max(pitchAngle, 5), 40);

    % --- Improved Bar Strength Calculation ---
    innerRegion = polarImg(1:ceil(end*0.25), :);
    innerFFT = fft2(innerRegion);
    innerFFTShift = fftshift(innerFFT);
    centerRowInner = ceil(size(innerFFTShift, 1)/2);
    centerColInner = ceil(size(innerFFTShift, 2)/2);

    bandWidth = 5;
    rowRange = max(centerRowInner - floor(bandWidth/2), 1):min(centerRowInner + floor(bandWidth/2), size(innerFFTShift, 1));
    innerProfile = mean(abs(innerFFTShift(rowRange, :)), 1);

    m0Component = innerProfile(centerColInner);

    m2IdxLeft = centerColInner - 2;
    m2IdxRight = centerColInner + 2;

    if m2IdxLeft > 0 && m2IdxRight <= length(innerProfile)
        m2Component = (innerProfile(m2IdxLeft) + innerProfile(m2IdxRight)) / 2;
    else
        m2Component = 0;
    end

    if m0Component > 0
        barStrength = m2Component / (m0Component + m2Component);
    else
        barStrength = 0;
    end
    barStrength = min(barStrength * 2, 1);
end

analyzeGalaxySpiralManualCut
