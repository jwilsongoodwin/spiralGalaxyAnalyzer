%% Real Example
function analyzeGalaxySpiral()
    % Open file selection dialog to choose a galaxy image
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.tif;*.bmp', 'Image Files (*.jpg, *.png, *.tif, *.bmp)'}, 'Select a Galaxy Image');
    
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('User cancelled file selection');
        return;
    end
    
    % Load the image
    fullPath = fullfile(pathname, filename);
    originalImg = imread(fullPath);
    
    % Display status
    disp('Image loaded. Processing...');
    
    % Convert to grayscale if image is RGB
    if size(originalImg, 3) > 1
        grayImg = rgb2gray(originalImg);
    else
        grayImg = originalImg;
    end
    
    % Convert to double for processing
    grayImg = im2double(grayImg);
    
    % Resize large images for faster processing
    if max(size(grayImg)) > 1200
        grayImg = imresize(grayImg, 1200/max(size(grayImg)));
        disp('Image resized to improve performance.');
    end
    
    % Apply logarithmic stretch to enhance spiral structure
    % This is crucial for better detection of spiral features
    grayImg = log1p(grayImg);
    grayImg = (grayImg - min(grayImg(:))) / (max(grayImg(:)) - min(grayImg(:)));
    
    % Display original image
    figure('Name', 'Galaxy Analysis', 'Position', [100, 100, 1200, 800]);
    subplot(2, 3, 1);
    imshow(originalImg);
    title('Original Galaxy Image');
    drawnow;
    
    % Display status
    disp('Finding galaxy center...');
    
    % Center detection and galaxy centering
    [centerX, centerY] = findGalaxyCenter(grayImg);
    
    % Allow manual center adjustment if needed
    choice = questdlg('Use detected center or manually specify?', 'Center Selection', ...
                      'Use Detected', 'Specify Manually', 'Use Detected');
    
    if strcmp(choice, 'Specify Manually')
        disp('Click on the galaxy center:');
        figure;
        imshow(grayImg);
        [centerX, centerY] = ginput(1);
        close;
    end
    
    % Display centered image with center marked
    subplot(2, 3, 2);
    imshow(grayImg);
    hold on;
    plot(centerX, centerY, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    title('Galaxy Center');
    drawnow;
    
    % Crop or pad to ensure square dimensions for polar transform
    [rows, cols] = size(grayImg);
    maxDim = max(rows, cols);
    squareImg = zeros(maxDim, maxDim);
    
    % Place original image in center of square canvas
    rowOffset = floor((maxDim - rows) / 2);
    colOffset = floor((maxDim - cols) / 2);
    squareImg((1:rows) + rowOffset, (1:cols) + colOffset) = grayImg;
    
    % Recalculate center coordinates for the square image
    centerX = centerX + colOffset;
    centerY = centerY + rowOffset;
    
    % Calculate maximum radius to use in polar transform
    maxRadius = min([centerX, centerY, maxDim-centerX, maxDim-centerY]);
    
    % Display status
    disp('Converting to polar coordinates...');
    
    % Create polar transform of the galaxy - HIGHER RESOLUTION FOR ACCURACY
    [polarImg, rVals, thetaVals] = cartesian2PolarImageFast(squareImg, centerX, centerY, maxRadius);
    
    % Display polar transform
    subplot(2, 3, 3);
    imshow(polarImg, []);
    title('Polar Transform (r, \theta)');
    xlabel('\theta (radians)');
    ylabel('Radius (pixels)');
    drawnow;
    
    % Display status
    disp('Calculating FFT...');
    
    % Apply 2D FFT to the polar image
    fftPolar = fft2(polarImg);
    fftShift = fftshift(fftPolar);
    logPower = log(1 + abs(fftShift));
    
    % Display FFT magnitude
    subplot(2, 3, 4);
    imshow(logPower, []);
    title('2D FFT of Polar Transform');
    colormap jet;
    drawnow;
    
    % Display status
    disp('Analyzing spiral parameters...');
    
    % Analyze the FFT to extract spiral parameters using improved methods
    [pitchAngle, barStrength, dominantMode] = analyzeImprovedSpiralFFT(fftShift, polarImg);
    
    % Display radial profile for pitch angle analysis with detected features
    subplot(2, 3, 5);
    [rows, cols] = size(fftShift);
    centerRow = ceil(rows/2);
    profile = abs(fftShift(centerRow, :));
    
    % Plot profile with detected peaks
    [pks, locs] = findpeaks(profile, 'SortStr', 'descend', 'NPeaks', 8);
    plot(profile, 'LineWidth', 1.5);
    hold on;
    
    % Highlight the dominant mode
    centerCol = ceil(cols/2);
    dominantPeakIdx = find(abs(locs - centerCol) == dominantMode, 1);
    
    if ~isempty(dominantPeakIdx)
        dominantLoc = locs(dominantPeakIdx);
        plot(dominantLoc, profile(dominantLoc), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        text(dominantLoc, profile(dominantLoc)*1.1, ['m = ', num2str(dominantMode)], 'FontSize', 12);
    end
    
    % Plot other peaks
    plot(locs, pks, 'go', 'MarkerSize', 6);
    
    title('Spectral Profile - Dominant Mode Detection');
    xlabel('Frequency Component');
    ylabel('Magnitude');
    grid on;
    
    % Display results
    subplot(2, 3, 6);
    axis off;
    text(0.1, 0.8, ['Estimated Pitch Angle: ' num2str(pitchAngle, '%.2f') ' degrees'], 'FontSize', 14);
    text(0.1, 0.7, ['Dominant Mode (arms): ' num2str(dominantMode)], 'FontSize', 14);
    text(0.1, 0.6, ['Bar Strength: ' num2str(barStrength, '%.2f')], 'FontSize', 14);
    
    % Add more detailed analysis
    if pitchAngle < 10
        spiralType = 'Tightly Wound (Sa-type)';
    elseif pitchAngle < 20
        spiralType = 'Intermediate (Sb-type)';
    else
        spiralType = 'Loosely Wound (Sc-type)';
    end
    
    text(0.1, 0.5, ['Galaxy Classification: ' spiralType], 'FontSize', 14);
    text(0.1, 0.3, ['Image: ' filename], 'FontSize', 12);
    title('Analysis Results');
    
    % Add overall title
    sgtitle(['Spiral Galaxy Analysis: ' filename], 'FontSize', 16);
    
    % Save results to a text file
    resultFile = fullfile(pathname, [filename(1:end-4), '_analysis.txt']);
    fid = fopen(resultFile, 'w');
    fprintf(fid, 'Spiral Galaxy Analysis Results\n');
    fprintf(fid, '----------------------------\n');
    fprintf(fid, 'Image: %s\n', filename);
    fprintf(fid, 'Pitch Angle: %.2f degrees\n', pitchAngle);
    fprintf(fid, 'Dominant Mode (arms): %d\n', dominantMode);
    fprintf(fid, 'Bar Strength: %.2f\n', barStrength);
    fprintf(fid, 'Galaxy Classification: %s\n', spiralType);
    fprintf(fid, '----------------------------\n');
    fprintf(fid, 'Analysis performed on: %s\n', datestr(now));
    fclose(fid);
    
    disp(['Analysis complete! Results saved to: ', resultFile]);
end

function [centerX, centerY] = findGalaxyCenter(img)
    % Improved center detection using both intensity and symmetry
    
    % Apply Gaussian smoothing to reduce noise
    smoothImg = imgaussfilt(img, 5);
    
    % Find the intensity-weighted centroid
    [rows, cols] = size(smoothImg);
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Find brightest region (likely the bulge)
    % Use top 5% of pixels for better accuracy
    threshold = quantile(smoothImg(:), 0.95);
    mask = smoothImg > threshold;
    
    if sum(mask(:)) > 0
        centerX = sum(X(:) .* mask(:) .* smoothImg(:)) / sum(mask(:) .* smoothImg(:));
        centerY = sum(Y(:) .* mask(:) .* smoothImg(:)) / sum(mask(:) .* smoothImg(:));
    else
        % Fall back to standard centroid if threshold method fails
        totalWeight = sum(smoothImg(:));
        centerX = sum(X(:) .* smoothImg(:)) / totalWeight;
        centerY = sum(Y(:) .* smoothImg(:)) / totalWeight;
    end
    
    % Round to nearest pixel
    centerX = round(centerX);
    centerY = round(centerY);
end

function [polarImg, rVals, thetaVals] = cartesian2PolarImageFast(img, centerX, centerY, maxRadius)
    % Convert a Cartesian image to polar coordinates - OPTIMIZED VERSION
    % Higher resolution for better spiral analysis
    
    % Define resolution for polar conversion
    numRadialBins = 200;  % Higher radial resolution
    numAngularBins = 360; % Full angular resolution (1 degree per bin)
    
    % Create polar coordinate arrays
    rVals = linspace(1, maxRadius, numRadialBins);
    thetaVals = linspace(0, 2*pi, numAngularBins+1);
    thetaVals = thetaVals(1:end-1); % Remove duplicate endpoint
    
    % Initialize polar image
    polarImg = zeros(numRadialBins, numAngularBins);
    
    % Create meshgrid for the entire transformation at once
    [THETA, R] = meshgrid(thetaVals, rVals);
    
    % Calculate corresponding Cartesian coordinates
    X = centerX + R .* cos(THETA);
    Y = centerY + R .* sin(THETA);
    
    % Use interpolation to get pixel values all at once
    polarImg = interp2(img, X, Y, 'linear', 0);
    
    % Apply logarithmic binning in radial direction to better visualize spiral pattern
    % This means more resolution near the center, less at the edges
    logRadialBins = exp(linspace(log(1), log(maxRadius), numRadialBins));
    [THETA_LOG, R_LOG] = meshgrid(thetaVals, logRadialBins);
    
    % Calculate corresponding Cartesian coordinates
    X_LOG = centerX + R_LOG .* cos(THETA_LOG);
    Y_LOG = centerY + R_LOG .* sin(THETA_LOG);
    
    % Create log-polar image (better for spiral analysis)
    logPolarImg = interp2(img, X_LOG, Y_LOG, 'linear', 0);
    
    % Combine regular and log-polar for better analysis
    polarImg = 0.5 * polarImg + 0.5 * logPolarImg;
end

function [pitchAngle, barStrength, dominantMode] = analyzeImprovedSpiralFFT(fftData, polarImg)
    % Analyze the FFT data to extract spiral parameters using improved methods

    % --- Initial Setup ---
    [rows, cols] = size(fftData);
    centerRow = ceil(rows/2);
    centerCol = ceil(cols/2);

    % --- Dominant Mode Detection (Improved) ---
    angularProfile = abs(fftData(centerRow, :));

    % Power spectrum around the center, excluding DC
    mValues = 1:10; % Try m = 1 to 10
    powerByMode = zeros(size(mValues));
    for i = 1:length(mValues)
        offset = mValues(i);
        leftIdx = centerCol - offset;
        rightIdx = centerCol + offset;
        if leftIdx > 0 && rightIdx <= cols
            powerByMode(i) = angularProfile(leftIdx) + angularProfile(rightIdx);
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

    % --- Spatial Domain: Polar Image (Log-polar Hough Transform) ---
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

    % Final weighted average
    pitchAngle = 0.7 * pitchAngle1 + 0.3 * pitchAngle2;
    pitchAngle = min(max(pitchAngle, 5), 40);

    % % M51 calibration correction
    % if dominantMode == 2 && pitchAngle > 30
    %     pitchAngle = pitchAngle * 0.6;
    % end

    % --- Improved Bar Strength Calculation ---
    innerRegion = polarImg(1:ceil(end*0.25), :);
    innerFFT = fft2(innerRegion);
    innerFFTShift = fftshift(innerFFT);
    centerRowInner = ceil(size(innerFFTShift, 1)/2);
    centerColInner = ceil(size(innerFFTShift, 2)/2);

    % Average across several rows for stability
    bandWidth = 5;
    rowRange = max(centerRowInner - floor(bandWidth/2), 1):min(centerRowInner + floor(bandWidth/2), size(innerFFTShift, 1));
    innerProfile = mean(abs(innerFFTShift(rowRange, :)), 1);

    m0Component = innerProfile(centerColInner);

    % Better m=2 detection with symmetric averaging
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
    barStrength = min(barStrength * 2, 1);  % Scaled to 0–1 range

end
analyzeGalaxySpiral