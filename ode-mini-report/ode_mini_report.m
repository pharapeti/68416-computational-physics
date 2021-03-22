clear; clc; close all

% Setup parameters
timestep = 0.01; % timestep (seconds)
totalTime = 100; % total time of simulation (seconds)
k = 0.25; % frequency term
gamma = 0.2; % dampening term

% Initial conditions
x_initial = 2; % initial position (metres)
v_initial = 0; % intial velocity (metres / second)
x0 = [x_initial; v_initial];

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
[t, x] = ode45(@(t, x) A * x, 0:timestep:totalTime, x0);
position = x(:, 1);
velocity = x(:, 2);

subplot(1, 3, 1);
plot(t, position);
title('Position vs Time');
xlabel('Time');
ylabel('Position');

subplot(1, 3, 2);
plot(t, velocity, 'r');
title('Velocity vs Time');
xlabel('Time');
ylabel('Velocity');

subplot(1, 3, 3);
plot(position, velocity, 'g', 'LineWidth', 1.2);
title('Postion vs Velocity');
xlabel('Position');
ylabel('Velocity');