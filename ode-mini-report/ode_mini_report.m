clear; clc; close all

% Setup parameters
timestep = 0.01; % timestep (seconds)
totalTime = 100; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;
k = 0.25; % frequency term
gamma = 0.25; % dampening term

% Initial conditions
x_initial = 2; % initial position (metres)
v_initial = 1; % intial velocity (metres / second)
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
%   Therefore -b^2 - gamma*b - k must equal 0
%   Solving for b
%   b = (gamma ± sqrt((-gamma)^2 - (4 * -1 * -k)) / -2
b1 = (gamma + sqrt((gamma)^2 - (4 .* k))) ./ -2;
b2 = (gamma - sqrt((gamma)^2 - (4 .*k))) ./ -2;

% Plot analytical solution
% figure(1);
% plot(timeSeries, analyticalSolution);

%% Part Two : Numberical Solution

% Convert 2nd order equation into a system of two coupled 1st order 
% equations using the variable suspension technique
%   Let xdot = v
%   Let vdot = xddot = -gamma * xdot - k * x
%
%   Since xdot = v, then...
%   d/dt [x; v] = [0 1; -k -gamma] * [x; v]

% Negatively dampened (will converge at y = 0)
A = [0 1; -k -gamma];

% Use ode45 to numerically solve
[t, x] = ode45(@(t, x) A * x, timeSeries, x0);
position = x(:, 1);
velocity = x(:, 2);

% Plot Position vs Time
figure(2);
subplot(1, 3, 1);
plot(t, position);
title('Position vs Time');
xlabel('Time');
ylabel('Position');

% Plot Velocity vs Time
subplot(1, 3, 2);
plot(t, velocity, 'r');
title('Velocity vs Time');
xlabel('Time');
ylabel('Velocity');

% Plot Position vs Velocity
subplot(1, 3, 3);
plot(position, velocity, 'g', 'LineWidth', 1.2);
title('Postion vs Velocity');
xlabel('Position');
ylabel('Velocity');