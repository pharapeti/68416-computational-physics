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

yFitted = modelFunction(p);

figure(1);
plot(x, y, '.', x, yFitted, '-');
title('Non-linear regression');
subtitle({ "where k = " + kOptimum, "and w = " + wOptimum});
legend('y', sprintf('%0.4f*exp(-%0.4f*x)', kOptimum, wOptimum));
xlabel('x');
ylabel('y');

% Generate domain for k and w values
k = linspace(-10, 10, 100);
w = linspace(-10, 10, 100).';
[kGrid, wGrid] = meshgrid(k, w);

% Generate meshgrid or broadcast to generate map of k and w values
grid = kGrid .* wGrid;

% Use map to compute the error function
% Store values in a 3 dimensional array, where dimenions are k, w, error
errorArray = nan([length(k), length(w)]);

% Not producing the correct result, we should expect the minimum point
% to be at kOptimum and wOptimum
for i = 1:length(k)
    for j = 1:length(w)
        errorArray(i, j) = errorFunction([k(i), w(j)]);
    end
end

% Use 2D map / contour to visualise the amplitude of error
figure(2);
imagesc(k, w, errorArray);

% Decorate Figure #1
colorbar
axis('image')
title('Error between model and dataset')
% subtitle(sprintf('%0.4f*exp(-%0.4f*x)', kOptimum, wOptimum));
xlabel('k')
ylabel('w')

% Highlight k/w values which minimise the errorFunction
lowestError = errorFunction([kOptimum, wOptimum]);

% Add lowestError to the 2D plot
