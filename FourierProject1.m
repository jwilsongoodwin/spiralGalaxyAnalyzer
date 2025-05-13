clear all
close all
%% 1D
x1 = 0:1:250;
t1 = 50;
t2 = 20;

Fs = 20;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
y1 = exp(-t./t1).*sin(2*2*pi*t);
y2 = exp(-t./t2).*sin(5*2*pi*t);
y3 = y1+y2;
X = y3 + randn(size(t));

figure()
plot(t, y1)
hold on
figure()
plot(t,y2)
figure()
plot(t,y3)
figure()
plot(t,X)
Y = fft(X);
figure()
plot(Fs/L*(0:L-1),abs(Y))
f = Fs/L*(0:(L/2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Z = ifft(Y)


% Plot the magnitude spectrum
figure(); 
plot(f, P1); 
title('Magnitude Spectrum of y3');
xlabel('Normalized Frequency'); 
ylabel('|Y(f)|');
grid on
figure()
plot(t,Z)

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
