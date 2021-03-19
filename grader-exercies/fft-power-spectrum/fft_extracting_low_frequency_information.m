clc; clear; close all;

t = linspace(0,1,1+50);
t = t(1:(end-1));
N = length(t);

% fit = c + bsin(wt) + acos(wt)
y = randn(1)+randn(1)*sin(2*pi*t)+randn(1)*cos(2*pi*t)+randn(size(t))*1e-2;

Y = fft(y); % do fft and extract components (remember to normalize)
Y = fftshift(y);
Y = Y / N;

% Reference: https://www.youtube.com/watch?v=qqqMK3fbLpU
% Determining coefficients
c = Y(1); % offset
b = -2 * imag(Y(2:N)); % imaginary component (sin)
a = -2 * real(Y(2:N)); % real component (cos)

% visually check your solution
plot(t, y, '.');

hold on;

fit = c + b .* sin(2*pi*t) + a .* cos(2*pi*t);
plot(t, fit, '-');