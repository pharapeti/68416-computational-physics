clc;
clear;
close all;

x = linspace(-1, 1, 1e3) * 3;
y = cos(x * 5 * (1+rand(1))) .* exp(-(x * (1+rand(1))) .^ 2) + randn(size(x)) * 0.05;

k = 0.12; % needs to be calculated from fitting yFit to y
w = 0.4;  % needs to be calculated from fitting yFit to y
yFit = cos(k * x) .* exp(-((w * x) .^2));

plot(x, y, '.', x, y, '-');
title('Function Fitting')
subtitle({ "where k = " + k, "and w = " + w});

latex = '\cos(kx)\exp(-(wx)^2)';
legend('Original Function', latex); % not rendering the latex

xlabel('x');
ylabel('y');