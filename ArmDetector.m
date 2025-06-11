function analyzeGalaxySpiralManualCut()
    % Select a galaxy image
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.tif;*.bmp', 'Images (*.jpg, *.png, *.tif, *.bmp)'}, 'Select a Galaxy Image');
    if isequal(filename, 0)
        disp('User cancelled.');
        return;
    end
    fullPath = fullfile(pathname, filename);
    imgRGB = imread(fullPath);

    % Convert to grayscale if needed
    if size(imgRGB, 3) > 1
        grayImg = rgb2gray(imgRGB);
    else
        grayImg = imgRGB;
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
    [polarImg, rVals, thetaVals] = cartesian2Polar(squareImg, centerX, centerY, maxRadius);

    % Normalize polar image before FFT
    polarImg = (polarImg - min(polarImg(:))) / (max(polarImg(:)) - min(polarImg(:)));

    % Preview polar image before FFT
    figure('Name', 'Polar Transformed Image');
    imshow(polarImg, []);
    title('Polar Image (Before FFT)');
    pause(2);

    % 2D FFT of polar image
    fftPolar = fft2(polarImg);
    fftShift = fftshift(fftPolar);

    % ----- DISPLAY 2D FFT magnitude and mark peak -----
    fftMag = abs(fftShift);
    [peakVal, peakIdx] = max(fftMag(:));
    [peakRow, peakCol] = ind2sub(size(fftMag), peakIdx);
    figure('Name', '2D FFT Magnitude (log scale)');
    imagesc(log(fftMag + 1));            % log-scale for better contrast
    axis image off;
    colormap(gray);
    hold on;
    plot(peakCol, peakRow, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    title(sprintf('FFT Magnitude (log)  –  Peak at (row=%d, col=%d)', peakRow, peakCol));
    pause(2);
 
    % ----------------------------------------------------

    % High-pass filter in Fourier space
    [rowsF, colsF] = size(fftShift);
    [X, Y] = meshgrid(1:colsF, 1:rowsF);
    centerRow = ceil(rowsF/2);
    centerCol = ceil(colsF/2);
    hpRadius = 3; % Reduced high-pass radius
    hpMask = sqrt((X - centerCol).^2 + (Y - centerRow).^2) > hpRadius;
    fftFiltered = fftShift .* hpMask;  % We'll hand this into the analysis routine

    % Analyze spiral (returns dominantMode, pitchAngle, barStrength)
    [pitchAngle, barStrength, dominantMode] = analyzeImprovedSpiralFFT(fftFiltered, polarImg);

    % ----- NEW: reconstruct and show which pixels belong to each arm -----
    showArmsOverlay( squareImg, imgRGB, centerX, centerY, rVals, thetaVals, dominantMode, fftShift );
    % ------------------------------------------------------------

    % Display results
    fprintf('Dominant Spiral Mode (Number of Arms): %d\n', dominantMode);
    fprintf('Estimated Pitch Angle: %.2f degrees\n', pitchAngle);
    fprintf('Bar Strength: %.2f\n', barStrength);
end

%% ------------------------------------------------------------------------
function [polarImg, rVals, thetaVals] = cartesian2Polar(img, centerX, centerY, maxRadius)
    numRadialBins   = 200;
    numAngularBins  = 360;
    rVals     = linspace(1, maxRadius, numRadialBins);
    thetaVals = linspace(0, 2*pi, numAngularBins);

    [THETA, R] = meshgrid(thetaVals, rVals);
    X = centerX + R .* cos(THETA);
    Y = centerY + R .* sin(THETA);

    polarImg = interp2(img, X, Y, 'linear', 0);
end

%% ------------------------------------------------------------------------
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
        leftIdx  = centerCol - offset;
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
    midPoint = floor(length(radialPhaseProfile)/2);
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
    centerRowInner = ceil(size(innerFFTShift,1)/2);
    centerColInner = ceil(size(innerFFTShift,2)/2);

    bandWidth = 5;
    rowRange = max(centerRowInner - floor(bandWidth/2), 1) : ...
               min(centerRowInner + floor(bandWidth/2), size(innerFFTShift,1));
    innerProfile = mean(abs(innerFFTShift(rowRange, :)), 1);

    m0Component = innerProfile(centerColInner);
    m2IdxLeft  = centerColInner - 2;
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

%% ------------------------------------------------------------------------
function showArmsOverlay(squareImg, imgRGB, centerX, centerY, rVals, thetaVals, m, fftShift)
    % Purpose: isolate the ±m Fourier bins, inverse-FFT back to polar, remap to Cartesian,
    %          then overlay on the original image so you see exactly which pixels
    %          belong to each of the m arms.

    % 1) Build a mask that keeps only the two bins at [centerCol ± m] on the centerRow
    [rowsF, colsF] = size(fftShift);
    centerRow = ceil(rowsF/2);
    centerCol = ceil(colsF/2);

    % A mask of zeros, then we copy over only the ±m bins:
    fftMask = zeros(size(fftShift));
    posIdx = centerCol + m;
    negIdx = centerCol - m;
    if posIdx > colsF, posIdx = posIdx - colsF; end
    if negIdx < 1,   negIdx = negIdx + colsF; end

    fftMask(centerRow, posIdx) = fftShift(centerRow, posIdx);
    fftMask(centerRow, negIdx) = fftShift(centerRow, negIdx);

    % 2) Inverse-shift and inverse-FFT to get a “filtered polar” image (only those m-arm frequencies)
    invShift = ifftshift(fftMask);
    filteredPolar = real(ifft2(invShift));

    % 3) Now we need to map filteredPolar (size: [numRadialBins x numAngularBins]) back to Cartesian.
    %    Create an empty Cartesian canvas the same size as squareImg, then for each (r,theta) bin,
    %    place the filteredPolar value at the corresponding (x,y) location:
    [numRadialBins, numAngularBins] = size(filteredPolar);
    [THETA, R] = meshgrid(thetaVals, rVals);  % THETA = [numRadialBins x numAngularBins], same as cartesian2Polar
    Xr = centerX + R .* cos(THETA);
    Yr = centerY + R .* sin(THETA);

    armCanvas = zeros(size(squareImg));  % will hold the “reconstructed” arms in Cartesian

    % Round to nearest pixel. Some points may map outside integer pixel grid; skip those.
    for i = 1:numRadialBins
        for j = 1:numAngularBins
            xPix = round(Xr(i,j));
            yPix = round(Yr(i,j));
            if xPix >= 1 && xPix <= size(armCanvas,2) && yPix >= 1 && yPix <= size(armCanvas,1)
                armCanvas(yPix, xPix) = filteredPolar(i,j);
            end
        end
    end

    % 4) Normalize armCanvas & create a simple transparency mask
    armCanvas = armCanvas - min(armCanvas(:));
    if max(armCanvas(:)) > 0
        armCanvas = armCanvas / max(armCanvas(:));
    end

    % 5) Overlay on the original (RGB) image
    figure('Name', sprintf('Overlay of %d Spiral Arm(s)', m));
    imshow(imgRGB, []);
    hold on;

    % Create a colored contour for arm pixels (threshold at, say, 20% of max)
    thresh = 0.2;
    [contourY, contourX] = find(armCanvas > thresh);
    plot(contourX, contourY, '.', 'MarkerSize', 4, 'Color', [1 0 0]); 
    % Red dots over pixels that strongly belong to the m arm

    title(sprintf('Pixels Belonging to the %d-Arm Mode (red)', m));
    hold off;
    % Pause so user can inspect
    pause(3);
end
analyzeGalaxySpiralManualCut;