clc;
clear;
close all;

% Define given dataset
x = linspace(-1, 1, 1e3) * 3;
y = cos(x * 5 * (1+rand(1))) .* exp(-(x * (1+rand(1))) .^ 2) ...
    + randn(size(x)) * 0.05;

% Define a lambda for model function which behaves like the actual data
% ...where p(1) is k
% ...where p(2) is w
modelFunction = @(p) cos(p(1) * x) .* exp(-((p(2) * x) .^2));

% Generate lamda for error function which returns the error between the
% model function and the actual data (based on the parameters passed in)
errorFunction = @(p) norm(y - modelFunction(p));

% Minimise error function providing an appropriate parameter estimate
p = fminsearch(errorFunction, [1 2]);

% Upack parameters returned by fminsearch which define the parameters that
% minimise the error function
kOptimum = p(1);
wOptimum = p(2);

% Generate fitted function
yFitted = modelFunction(p);

% Plot and decorate the actual dataset, and the modelled function
figure(1);
plot(x, y, '.', x, yFitted, '-');
title('Non-linear regression');
legend('y', sprintf('%0.4f*exp(-%0.4f*x)', kOptimum, wOptimum));
xlabel('x');
ylabel('y');

%% Part Two
% Generate domain for k and w values
k = linspace(-10, 10, 100);
w = linspace(-10, 10, 100);

% Generate error over the k, w domain
errorArray = nan([length(k), length(w)]);
for i = 1:length(k)
    for j = 1:length(w)
        kVal = k(:, i);
        wVal = w(:, j);

        errorArray(i, j) = errorFunction([kVal, wVal]);
    end
end

% Use 2D map / contour to visualise the amplitude of error
figure(2);
imagesc(w, k, errorArray);
hold on;

% Highlight k/w value which minimises the error function
lowestError = errorFunction([kOptimum, wOptimum]);

% Add lowestError to the 2D plot
plot3(wOptimum, kOptimum, 0, '-s', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'red');

% Decorate Figure #1
colorbar
axis('image')
title('Error between model and dataset')
subtitle('Smaller the error, the more accurate the fit');
xlabel('w')
ylabel('k')
