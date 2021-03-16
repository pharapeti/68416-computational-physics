clear; clc; close all;

% Calculate and plot the absolute of fft for the specified signal.
% Is the peak frequency correct, and can you explain why the spectrum looks
% this way?

fs = 1e3; % Sampling frequency in Hertz
t = 0:(1/fs):1; % Sampling period
t = t(1:(end-1)); % Time vector
y = sin(2*pi*100*t); % Signal

subplot(1, 2, 1);
plot(t, y);
title('Signal');

Y = fft(y); % FFT computation
Y = fftshift(Y); % unfold

dt = mean(diff(t)); % sample spacing
N = length(t);
df = 1/(N*dt); % frequency spacing
fi = (0:(N-1)) - floor(N/2); % unfolded index
f = df*fi; % frequency vector

subplot(1, 2, 2);
plot(f, abs(Y/N)); % amplitude vs frequency
% plot(f, abs(Y).^2); % power vs frequency
title('FFT');