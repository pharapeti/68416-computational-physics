clc;
clear;
close all;

x = linspace(0,(1+rand(1)))*5;
y = (1+rand(1))*exp(-(1+rand(1))*x)+randn(size(x))*1e-2;

% put the model function here: p will contain A and k 
f = @(p) p(1) .* exp(-p(2) .* x);

% put the merit function here
% which returns a single value describing the 
% difference between model and data
m = @(p) norm(y - f(p));

% optimize: put in the merit function and starting guess for p
p = fminsearch(m, [1 2]);

% unpack
A = p(1);
k = p(2);

% fitted y goes here (use the model function)
yf= f(p);

plot(x, y, '.', x, yf, '-') % plot with overlay
legend('y', sprintf('%0.4f*exp(-%0.4f*x)', A, k));