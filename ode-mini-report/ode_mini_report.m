%% ODE Mini-Report Assignment (By Patrice Harapeti)

%% Background
% Damped Simple Oscillator

%% Part Zero : Setup

% Clear existing workspace
clear; clc; close all

% Setup parameters
timestep = 0.01; % timestep (seconds)
totalTime = 100; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;
k = 0.25; % frequency term
gamma = 0.25; % dampening term

% Initial conditions
x_initial = 2; % initial position (metres)
v_initial = 0; % intial velocity (metres / second)
x0 = [x_initial; v_initial];

%% Part One : Analytical Solution

%   Derivation
%   Let x = e^bt
%   Let xdot = b * e^bt
%   Let xddot = b^2 * e^bt
%
%   Therefore, -xddot - gamma*xdot - k*x = 0 equals
%   -b^2 * e^bt - (gamma * b * e^bt) - ke^bt = 0
%
%   Pull e^bt out as common factor, this leaves us
%   e^bt (-b^2 - gamma*b - k) = 0
%
%   Therefore b^2 + gamma*b + k must equal 0
%   Solving for b
%   b = (gamma ± sqrt((gamma)^2 - (4k)) / -2
%   
%   Overdamped when...         gamma^2 - 4k > 0
%   Critically damped when...  gamma^2 - 4k = 0
%   Underdamped when...        gamma^2 - 4k < 0

% Calculate roots of characteristics equations
b_1 = (-gamma + sqrt(gamma.^2 - (4 .* k))) ./ 2;
b_2 = (-gamma - sqrt(gamma.^2 - (4 .* k))) ./ 2;

% Solve for A and B constants
A = (-2 .* b_2) ./ (b_1 - b_2);
B = (2 .* b_1) ./ (b_1 - b_2);

alpha = real(b_1);
beta = imag(b_2);

% Define function in three cases based on the determinant of the roots
% Reference: https://nrich.maths.org/11054
over_damped = (A .* exp(b_1 .* timeSeries)) + (B .* exp(b_2 .* timeSeries));
critically_damped = (A + B.*timeSeries) .* exp(b_1 .* timeSeries);
under_damped = exp(alpha .* timeSeries) .* (A.*cos(beta.*timeSeries) + B.*sin(beta.*timeSeries));

% Plot analytical solution
figure(1);
subplot(1, 3, 1);
plot(timeSeries, over_damped);
title('Overdamped');
xlim([0, 70]);
xlabel('Time');
ylabel('Position');

subplot(1, 3, 2);
plot(timeSeries, critically_damped);
title('Critically Damped');
xlim([0, 70]);
xlabel('Time');
ylabel('Position');

subplot(1, 3, 3);
plot(timeSeries, under_damped);
title('Underdamped');
xlim([0, 70]);
xlabel('Time');
ylabel('Position');

%% Part Two : Numberical Solution (USE EULERS OR ANOTHER METHOD)

% Convert 2nd order equation into a system of two coupled 1st order 
% equations using the variable suspension technique
%   Let xdot = v
%   Let vdot = xddot = -gamma * xdot - k * x
%
%   Since xdot = v, then...
%   d/dt [x; v] = [0 1; -k -gamma] * [x; v]

% Negatively dampened (will converge at y = 0)
% Also known as an underdamped system
A = [0 1; -k -gamma];

% Use ode45 to numerically solve
[t, x] = ode45(@(t, x) A * x, timeSeries, x0);
position = x(:, 1);
velocity = x(:, 2);

%% Part Three : Plotting

% Plot Position vs Time
figure(2);
subplot(1, 3, 1);
plot(t, position);

% Draw horizontal line at y = 0 to represent convergence value
yline(0, '--');
grid on;
title('Position vs Time');
xlim([0, 70]);
xlabel('Time');
ylabel('Position');

% Plot Velocity vs Time
subplot(1, 3, 2);
plot(t, velocity, 'r');
grid on;
title('Velocity vs Time');
xlim([0, 70]);
xlabel('Time');
ylabel('Velocity');

% Plot Position vs Velocity
subplot(1, 3, 3);
plot(position, velocity, 'g', 'LineWidth', 1.2);
grid on;
title('Postion vs Velocity');
xlabel('Position');
ylabel('Velocity');

% Plot surface plot of gamma/k vs time to visualise how the 
% ratio (dampeing) of the parameters affect the function

%% Part Four : Convergence

% Pick particular case
% Prove the numerical solution converges at the same value the analytical
% solution conveges to

%% Part Five : Error between analytical and numerical solution

% Pick particular case
% Produce a plot of the error between the analytical and numberical
% solution as as the step size of increased/decreased
% Plot should be plot(stepsize, error)

%% Part Six : Fourier Analysis of numerical solution

% Do this for each case
% Compute a power spectrum graph of the numerical solution
% Plot (frequency, power_freqSpace);

% Determine the center frequency and FWHM
% Include these in the plot