clc;
clear;
close all;

x = linspace(0,(1+rand(1)))*5; 
y = exp(-(1+rand(1))*x)+randn(size(x))*1e-2;

% put the model function here
% k should be approx. 1.45 (by inspection)
f = @(k) exp(-k .* x);

% put the merit function here
% which returns a single value describing the 
% difference between model and data
m = @(k) norm(y - f(k));

% optimize: put in the merit function and starting guess for k
p = fminsearch(m, 1.5);

% fitted y goes here (use the model function)
yf = f(p);
k = p;

% plot with overlay
plot(x, y, '.', x, yf, '-');
legend('y',sprintf('exp(-%0.4f*x)', k));