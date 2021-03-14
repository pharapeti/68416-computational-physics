clc; clear; close all;

% Given dataset
x = linspace(-1,1,1e3).' * 5;
y = exp(-x.^2) .* polyval(randn([1,5]),x) + randn(size(x)) * 0.02;

% Build your system matrix here
A = [x x x x];

% Solve the system, and then unpack the individual values into a and b
c = A \ y;
c1 = c(1);
c2 = c(2);
c3 = c(3);
c4 = c(4);

% Calculate your fitted y here using c
yf = exp(-x.^2) .* 1

% % Sum squared error
% S = sum((y - yf) .^ 2);
% 
% % Mean squared error
% % ... m = number of data points = 100
% % ... n = number of fit functions / coefficients = 2
% % ... therefore m - n
% m = 100;
% n = 2;
% 
% meanSquaredError = sqrt(S / (m - n));
% 
% % Uncertainty
% uncertainty = meanSquaredError .* sqrt(inv(A * A.'));
% 
% % % Estimate the errors in c
% dc = uncertainty * c;

plot(x, y, '.', x, yf, '-');
title('Linear regression');
% legend('y', sprintf('%0.4f*exp(-%0.4f*x)', kOptimum, wOptimum));
xlabel('x');
ylabel('y');