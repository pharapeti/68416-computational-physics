clc; clear; close all;

%% Part 1

% Given dataset
x = linspace(-1,1,1e3).' * 5;
y = exp(-x.^2) .* polyval(randn([1,5]),x) + randn(size(x)) * 0.02;

% Build model
A = exp(-x.^2) .* [x.^0 x.^1 x.^2 x.^3 x.^4];

% Solve the system for the coefficients for each power of x in the model
c = A \ y;

% Calculate your fitted function here using c and the model
yFit = A * c;

% Sum squared error
S = sum((y - yFit) .^ 2);

% Mean squared error
% ... m = number of data points = length(x)
% ... n = number of coefficients = length(c)
% ... therefore m - n
m = length(x);
n = length(c);
meanSquaredError = sqrt(S / (m - n));
fprintf('Mean Squared Error = %f\n', meanSquaredError);

% Generate Covariance Matrix from diaglon of inverse conjugate matrix of A
covarianceMatrix = diag(inv(A.' * A));

% Uncertainty Matrix
uncertainty = meanSquaredError .* sqrt(covarianceMatrix);
disp("Uncertainty for each coefficient of x = [" + ...
    num2str(uncertainty.') + ']');

figure(1);
plot(x, y, '.', 'LineWidth', 1);
hold on;
plot(x, yFit, '-', 'LineWidth', 2);
title('Linear regression');
xlabel('x');
ylabel('y');
subtitle("Coefficients of x = [" + num2str(c.') + ']');
legend({'data', '$exp(-x^2)\sum\limits_{i=0}^4c_{i}x^{i}$'}, ...
    'Interpreter', 'Latex');

%% Part 2

% Delete temporarily variables of Part 2 (memory optimisation)
clear;

% Specify given dataset
x = linspace(-1,1,1+1e3).' * 10;
y = exp(-(x/(1+rand(1))).^2) .* cos(2*pi*(1+rand(1))*x);

% Use fft function to generate Discrete Fourier Transform (DFT)
Y = fft(y);

% Shift the zero-frequency component to the center of the array
Y = fftshift(Y);

% Normalise amplitude by normalisation factor, due to summation in DFT
Y = Y / sqrt(length(x));

% Calculate power spectrum of signal
power = abs(Y) .^ 2;

% Limit power to be positive for frequency analysis
xPositive = x(floor(length(x)/2) + 1:length(x));
powerPositive = power(floor(length(power)/2)+1:length(power));

% Plot Amplitude vs x
figure(2);
subplot(1, 2, 1);
plot(x, y);
title('Amplitude vs x');
xlim([-5, 5]);
xlabel('x');
ylabel('Amplitude');

% Plot Power vs Frequency
subplot(1, 2, 2);
plot(x, power, '-r');
title('Power vs Frequency');
xlim([-1, 1]);
xlabel('Hz');
ylabel('Power');

% Plot Power vs Frequency for positive frequencies
figure(3);
plot(xPositive, powerPositive, '-r');
title('Power vs Frequency');
subtitle('For positive frequencies');
xlabel('Hz');
ylabel('Power');
hold on;

% Expectation values
% E (hz) = Sum of p(x) * x / Sum of p(x)
standardDeviation = std(powerPositive);
% expectationValues = (sum(powerPositive) .* 2) ./ sum(powerPositive);

% Full-width half maximum of positive power
% FWHM = expectationValues .^2 / sqrt(log(2));
FWHM = 2 .* sqrt(2 .* log(2)) .* standardDeviation;

% Find maxima of positive power vs frequency
[powerMaxima, maximaIndex] = max(powerPositive);
frequencyMaxima = xPositive(maximaIndex);

% Mark maxima on plot of positive power vs frequency
plot(frequencyMaxima, powerMaxima, 'xb', 'LineWidth', 2);
hold on;

% Calculate half-max power threshold
powerThreshold = powerMaxima / 2;

% Find indices where power exceeds threshold
indices = find(powerPositive > powerThreshold);

% Calculate width at half-maximum
fullWidth = xPositive(indices(end)) - xPositive(indices(1));

% Calculate center of positive power vs frequency
centerIndex = xPositive(floor(length(indices)/2) + 1);

plot(xPositive(indices), powerThreshold * ones(size(xPositive(indices))),...
    '--', 'LineWidth', 2);

% Center plot around center value using arbitrary distance from center
% Based the offset on the full width of the half-maximum as I believe
% it's fair to assume the Gaussian distribution does not extend any further
% than 3 times the peak in either direction 
xlim([frequencyMaxima - (3 * fullWidth), ...
    frequencyMaxima + (3 * fullWidth)]);
