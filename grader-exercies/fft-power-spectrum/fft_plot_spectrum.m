clear; clc; close all;

% Calculate and plot the absolute of fft for the specified signal.
% Is the peak frequency correct, and can you explain why the spectrum looks
% this way?

fs = 1e3; % Sampling frequency in Hertz
t = 0:(1/fs):1; % Sampling period
t = t(1:(end-1)); % Time vector
N = length(t); % Number of samples
y = sin(2*pi*100*t); % Signal

subplot(1, 3, 1);
plot(t, y);
title('Signal');
xlabel('Time (s)');
ylabel('Amplitude');

Y = fft(y); % FFT computation
Y = fftshift(Y); % unfold
Y = Y/N;

dt = mean(diff(t)); % sample spacing
df = 1 / (N * dt); % frequency spacing
fi = (0:(N-1)) - floor(N/2); % unfolded index
f = df * fi; % frequency vector

subplot(1, 3, 2);
plot(f, abs(Y)); % amplitude vs frequency
title('FFT');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(1, 3, 3);
plot(f, abs(Y).^2); % power vs frequency
title('FFT');
xlabel('Frequency (Hz)');
ylabel('Power');