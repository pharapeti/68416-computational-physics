clc; clear; close all;

% Decay rate
% dN (t) / dT = - k * N(t)

% Step size
dT = 0.01;

% Time
time = 0:dT:100;

% Population array
population = nan(size(time));
population_inital = 10000;
population(1) = population_inital;

% Decay constant
k = 0.02;

% Find population at time t = 100s
for i = 1:length(population) - 1
    population(i + 1) = population(i) + (-k .* dT .* population(i));
end

plot(time, population);