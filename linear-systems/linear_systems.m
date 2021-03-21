clc; clear; close all;

%% Part 1

% Given dataset
x = linspace(-1,1,1e3).' * 5;
y = exp(-x.^2) .* polyval(randn([1,5]),x) + randn(size(x)) * 0.02;

% Build your system matrix here
A = [x.^0 x.^1 x.^2 x.^3 x.^4];

% Not sure if I need to do this!!!
% Rearrange to produce y/exp(-x.^2) = A * c
%yMod = y ./ exp(-x.^2);

% Solve the system, and then unpack the values of the coefficients for each
% power of x in the summation
c = A \ y;
c0 = c(1);
c1 = c(2);
c2 = c(3);
c3 = c(4);
c4 = c(5);

% Calculate your fitted y here using c
% yf = exp(-x.^2) .* y;

yFit = A * c;
yFit = yFit .* exp(-x.^2);

% Sum squared error
% S = sum((y - yFit) .^ 2);

% Mean squared error
% ... m = number of data points = 100
% ... n = number of fit functions / coefficients = 2
% ... therefore m - n
% m = 100;
% n = 2;
% meanSquaredError = sqrt(S / (m - n));

% Uncertainty
% uncertainty = meanSquaredError .* sqrt(inv(A * A.'));

% Estimate the errors in c
% dc = uncertainty * c;

%figure(1);
% plot(x, y, '.', x, yFit, '-');
%title('Linear regression');
%xlabel('x');
%ylabel('y');
% legend('Dataset', {'$exp(-x^2)\sum\limits_{i=0}^4c_{i}x^{i}$'}, 'Interpreter', 'Latex');

%% Part 2

% Given dataset
x = linspace(-1,1,1+1e3).' * 10;
y = exp(-(x/(1+rand(1))).^2) .* cos(2*pi*(1+rand(1))*x);

Y = fft(y);
Y = fftshift(Y);
Y = Y / length(x);

power = abs(Y) .^ 2;

figure(2);
subplot(1, 2, 1);
plot(x, y);
title('Amplitude vs Frequency');
xlabel('Hz');
ylabel('Amplitude');

subplot(1, 2, 2);
plot(x, power, '-r');
title('Power vs Frequency');
xlabel('Hz');
ylabel('Power');