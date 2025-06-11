function analyze_spiral_galaxy_logpolar()
% analyze_spiral_galaxy_logpolar
%   Load any spiral‐galaxy image, let user choose (or auto‐detect) the center,
%   mask out the core, convert to log‐polar coordinates (u = ln(r)),
%   compute a 2D FFT in (u, θ) space, and then estimate the number of arms (m)
%   and the pitch angle.
%
% USAGE:
%   1) Place this file in your MATLAB path.
%   2) At the MATLAB prompt, type: analyze_spiral_galaxy_logpolar
%   3) Follow the prompts in the figure windows.

    %% 1) Read in an image
    [fname, path] = uigetfile({'*.jpg;*.png;*.tif;*.bmp','Galaxy images (*.jpg, *.png, *.tif, *.bmp)'}, ...
                              'Select a spiral galaxy image');
    if isequal(fname,0)
        disp('No file selected. Exiting.');
        return;
    end
    img_rgb = imread(fullfile(path, fname));

    % Convert to grayscale if necessary
    if size(img_rgb,3) == 3
        img_gray = im2double(rgb2gray(img_rgb));
    else
        img_gray = im2double(img_rgb);
    end

    figure('Name','Original Image','NumberTitle','off');
    imshow(img_gray,[]);
    title('Original Galaxy Image');

    %% 2) Choose center: manual or automatic
    choice = questdlg('How would you like to select the galaxy center?', ...
                      'Center Selection', ...
                      'Manual (click)','Automatic (intensity centroid)','Manual (click)');
    switch choice
        case 'Manual (click)'
            title('Click on the galaxy center');
            [x0, y0] = ginput(1);
            x0 = round(x0); 
            y0 = round(y0);
        case 'Automatic (intensity centroid)'
            % Binarize and find centroid of the largest region
            bw = im2bw(imadjust(img_gray), graythresh(img_gray));
            props = regionprops(bw, 'Area', 'Centroid');
            if isempty(props)
                [ny, nx] = size(img_gray);
                x0 = round(nx/2); 
                y0 = round(ny/2);
                warning('Automatic thresholding failed; using image center.');
            else
                % Pick region with maximum area
                areas = [props.Area];
                [~, idx_max] = max(areas);
                centroid = props(idx_max).Centroid;
                x0 = round(centroid(1)); 
                y0 = round(centroid(2));
            end
        otherwise
            disp('No selection—exiting.');
            return;
    end
    hold on;
    plot(x0, y0, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;

    %% 3) Mask the core with a circular ROI
    promptMsg = 'Draw a circle around the galaxy core to mask it (double‐click to finish).';
    hFig = figure('Name','Mask the Core','NumberTitle','off');
    imshow(img_gray,[]);
    title(promptMsg);
    % Initial small circle around (x0,y0)
    initRadius = 20;
    hEllipse = imellipse(gca, [x0-initRadius, y0-initRadius, 2*initRadius, 2*initRadius]);
    position = wait(hEllipse);  % Wait until user double-clicks
    maskCore = createMask(hEllipse);
    close(hFig);

    % Invert mask: 1 outside core, 0 inside core
    maskOutside = ~maskCore;
    img_masked = img_gray .* maskOutside;

    figure('Name','Masked Image','NumberTitle','off');
    imshow(img_masked,[]);
    title('Galaxy with Core Masked');

    %% 4) Convert to log‐polar coordinates (u = ln(r))
    [ny, nx] = size(img_masked);
    % Choose r_min just outside the masked core
    % Estimate r_min as the smallest nonzero distance from (x0,y0) to masked region
    distMap = sqrt((repmat((1:nx), ny,1) - x0).^2 + (repmat((1:ny)',1,nx) - y0).^2);
    distMap(~maskOutside) = Inf;  % ignore core
    r_min = max(1, floor(min(distMap(distMap<Inf))));
    % r_max must fit within image bounds
    r_max = floor(min([x0, y0, nx - x0, ny - y0]));

    % Sampling parameters
    nu = 512;       % number of u samples (log‐radius)
    ntheta = 512;   % number of θ samples
    uVec = linspace(log(r_min), log(r_max), nu);
    thetaVec = linspace(0, 2*pi, ntheta);
    [U, Theta] = meshgrid(uVec, thetaVec);

    R = exp(U);
    Xp = R .* cos(Theta) + x0;
    Yp = R .* sin(Theta) + y0;
    img_logpolar = interp2(img_masked, Xp, Yp, 'linear', 0);

    figure('Name','Log‐Polar Mapped','NumberTitle','off');
    imagesc(uVec, thetaVec, img_logpolar);
    axis xy;
    xlabel('u = ln(r)');
    ylabel('θ (radians)');
    title('Log‐Polar (u, θ) Intensity');
    colormap gray;

    %% 5) Subtract mean along each θ to remove m=0 component
    img_lp_zeroMean = img_logpolar - mean(img_logpolar, 2);

    %% 6) Compute 2D FFT in (u, θ) space
    F_lp = fft2(img_lp_zeroMean);
    F_lp_shift = fftshift(F_lp);
    P_lp = abs(F_lp_shift);

    figure('Name','2D FFT (u, θ)','NumberTitle','off');
    imagesc((-ntheta/2):(ntheta/2-1), (-nu/2):(nu/2-1), log10(P_lp + eps));
    axis xy;
    xlabel('Angular frequency m');
    ylabel('Radial frequency k_u');
    title('Log‐Power of 2D FFT in (u, θ)');
    colormap jet;
    colorbar;

    %% 7) Find dominant m_peak and k_u_peak
    half_m = floor(ntheta/2);
    half_k = floor(nu/2);
    m_vals = 1:(half_m-1);

    power_vs_m = zeros(size(m_vals));
    for idx = 1:length(m_vals)
        col_index = half_m + m_vals(idx);
        power_vs_m(idx) = sum(P_lp(:, col_index));
    end
    [~, imax] = max(power_vs_m);
    m_peak = m_vals(imax);

    col_index_peak = half_m + m_peak;
    rows_to_search = (half_k+2):(nu-1);  % skip k_u = 0 at row half_k+1
    [~, rel_idx] = max(P_lp(rows_to_search, col_index_peak));
    ir_peak = rows_to_search(rel_idx);
    k_u_peak = abs(ir_peak - (half_k+1));

    % Compute pitch angle
    pitch_rad = atan(m_peak / k_u_peak);
    pitch_deg = pitch_rad * (180/pi);

    %% 8) Display results in Command Window
    fprintf('Estimated number of arms (|m|): %d\n', m_peak);
    fprintf('Estimated radial frequency in log‐r (k_u_peak): %d\n', k_u_peak);
    fprintf('Estimated pitch angle: %.2f°\n', pitch_deg);

    %% 9) Optional: Reconstruct single (m_peak, k_u_peak) mode and show
    mask_single = zeros(size(F_lp));
    F_unshift = ifftshift(F_lp_shift);

    idx_m_pos = m_peak + 1;            % MATLAB indexing for m_peak > 0
    idx_k_pos = ir_peak;               % radial index in FFT array
    mask_single(idx_k_pos, idx_m_pos) = F_unshift(idx_k_pos, idx_m_pos);

    % Also include its complex‐conjugate counterpart (negative m)
    idx_m_neg = mod(ntheta - m_peak + 1, ntheta);
    if idx_m_neg == 0
        idx_m_neg = ntheta;
    end
    mask_single(idx_k_pos, idx_m_neg) = F_unshift(idx_k_pos, idx_m_neg);

    % Inverse FFT to get back to (u, θ)
    F2_single = fftshift(mask_single);
    img_mode = real(ifft2(F2_single));

    figure('Name','Dominant Mode Reconstruction','NumberTitle','off');
    imagesc(uVec, thetaVec, img_mode);
    axis xy;
    xlabel('u = ln(r)');
    ylabel('θ (radians)');
    title(sprintf('Reconstructed (m = %d) Mode in (u, θ)', m_peak));
    colormap parula;
    colorbar;

    disp(' ');
    disp('Done. Close figures when finished viewing.');
end
