
clear all
close all
clc

%% ----------------- 1D Fourier Transform -----------------
x1 = 0:1:250;
t1 = 50;
t2 = 20;

Fs = 20;           % Sampling frequency
T = 1/Fs;          % Sampling period
L = 1000;          % Length of signal
t = (0:L-1)*T;     % Time vector

% Signals
y1 = exp(-t./t1).*sin(2*2*pi*t);    % 2 Hz damped sine
y2 = exp(-t./t2).*sin(5*2*pi*t);    % 5 Hz damped sine
y3 = y1 + y2;                       % Combined signal
X = y3 + randn(size(t));           % Add Gaussian noise

% Plot y1
figure()
plot(t, y1)
xlabel('Time (s)')
ylabel('Amplitude')
title('y_1(t) = e^{-t/t1} \cdot sin(4\pi t)  [2 Hz]')
grid on

% Plot y2
figure()
plot(t, y2)
xlabel('Time (s)')
ylabel('Amplitude')
title('y_2(t) = e^{-t/t2} \cdot sin(10\pi t)  [5 Hz]')
grid on

% Plot combined signal y3
figure()
plot(t, y3)
xlabel('Time (s)')
ylabel('Amplitude')
title('Combined Signal y_3(t) = y_1 + y_2')
grid on

% Plot noisy signal
figure()
plot(t, X)
xlabel('Time (s)')
ylabel('Amplitude')
title('Noisy Signal X(t) = y_3 + noise')
grid on

% FFT
Y = fft(X);
f = Fs*(0:(L/2))/L;           % One-sided frequency axis
P2 = abs(Y/L);                % Two-sided spectrum
P1 = P2(1:L/2+1);             % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Adjust amplitude

% Plot magnitude spectrum and annotate frequency peaks
figure(); 
plot(f, P1); 
xlabel('Frequency (Hz)'); 
ylabel('|Y(f)|');
title('Magnitude Spectrum with Peak Annotations');
grid on
hold on

% Peak detection with threshold to avoid noise peaks
[peaks, locs] = findpeaks(P1, f, 'MinPeakHeight', max(P1)*0.3);

% Plot and label each detected peak
for i = 1:length(locs)
    % Draw a red circle at the peak
    plot(locs(i), peaks(i), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
    
    % Label the peak
    if abs(locs(i) - 2) < 0.5
        label = '2 Hz (y_1)';
    elseif abs(locs(i) - 5) < 0.5
        label = '5 Hz (y_2)';
    else
        label = sprintf('%.2f Hz', locs(i));
    end
    
    text(locs(i)+0.2, peaks(i), label, 'Color', 'red', ...
        'FontWeight', 'bold', 'FontSize', 10);
end


% Inverse FFT reconstruction
Z = ifft(Y);
figure()
plot(t, real(Z))
xlabel('Time (s)')
ylabel('Amplitude')
title('Reconstructed Signal from Inverse FFT of X(t)')
grid on


%% 2D
% Create checkerboard image
n = 8;  % number of checkers per row/column
squareSize = 32; % pixels per checker
imgSize = n * squareSize;

% Create checkerboard pattern
checker = checkerboard(squareSize, n, n) > 0.5;
checker = double(checker);  % convert to double precision image

% Show original checkerboard
figure;
imshow(checker);
title('Original Checkerboard Image');
% Compute the 2D FFT
F = fft2(checker);
F_shifted = fftshift(F);  % Center zero frequency

% Compute magnitude and phase spectra
magF = abs(F_shifted);
logMagF = log(1 + magF);  % for visualization

% Show log-magnitude of Fourier spectrum
figure;
imshow(logMagF, []);
title('Log Magnitude Spectrum (Centered)');

% Frequency axes (normalized from -0.5 to 0.5)
fx = linspace(-0.5, 0.5, imgSize);
fy = linspace(-0.5, 0.5, imgSize);

% Show frequency plot with axes
figure;
imagesc(fx, fy, logMagF);
axis image;
xlabel('fx'); ylabel('fy');
title('Calibrated Frequency Domain (Log Magnitude)');
colorbar;

% Inverse transform
reconstructed = ifft2(F);

% Visualize real part (should match original image)
figure;
imshow(real(reconstructed), []);
title('Reconstructed Image from Inverse FFT');

% Add a white circle to checkerboard
[x, y] = meshgrid(1:size(checker,2), 1:size(checker,1));  % Match image size
circle = sqrt((x - size(checker,2)/2).^2 + (y - size(checker,1)/2).^2) < 30;

% Add the circle to the checkerboard
checker_mod = checker + double(circle);  % Convert logical to double
% Show modified image
figure;
imshow(checker_mod, []);
title('Modified Image with Circle');

% Show modified Fourier transform
F_mod = fftshift(fft2(checker_mod));
figure;
imshow(log(1 + abs(F_mod)), []);
title('Fourier Transform of Modified Image');
% Add Gaussian noise
checker_noise = checker + 0.3*randn(size(checker));

% Show noisy image
figure;
imshow(checker_noise, []);
title('Noisy Checkerboard');

% Show FFT of noisy image
F_noise = fftshift(fft2(checker_noise));
figure;
imshow(log(1 + abs(F_noise)), []);
title('Fourier Transform of Noisy Image');
% Create low-pass filter (centered circle)
[rows, cols] = size(F_shifted);
[xx, yy] = meshgrid(1:cols, 1:rows);
mask_radius = 40;
low_pass_mask = sqrt((xx - cols/2).^2 + (yy - rows/2).^2) < mask_radius;

% Apply filter
F_filtered = F_shifted .* low_pass_mask;

% Inverse FFT
filtered_img = ifft2(ifftshift(F_filtered));

% Show result
figure;
imshow(real(filtered_img), []);
title('Low-Pass Filtered Image');
% Create a high-pass mask (inverse of low-pass)
[rows, cols] = size(F_shifted);
[xx, yy] = meshgrid(1:cols, 1:rows);
mask_radius = 40;
low_pass_mask = sqrt((xx - cols/2).^2 + (yy - rows/2).^2) < mask_radius;

% Invert it for high-pass
high_pass_mask = ~low_pass_mask;

% Apply the high-pass filter
F_highpass = F_shifted .* high_pass_mask;

% Inverse FFT
F_highpass_unshift = ifftshift(F_highpass);
img_highpass = real(ifft2(F_highpass_unshift));

% Display result
figure;
imshow(img_highpass, []);
title('High-Pass Filtered Image');
% FFT of low-pass filtered image
F_after_lowpass = fftshift(fft2(real(filtered_img)));
figure;
imshow(log(1 + abs(F_after_lowpass)), []);
title('Fourier Transform After Low-Pass Filtering');
% FFT of high-pass filtered image
F_after_highpass = fftshift(fft2(img_highpass));
figure;
imshow(log(1 + abs(F_after_highpass)), []);
title('Fourier Transform After High-Pass Filtering');
