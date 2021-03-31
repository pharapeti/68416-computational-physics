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

% Initial conditions
x_initial = 2; % initial position (metres)
v_initial = 0; % intial velocity (metres / second)
x0 = [x_initial; v_initial];

%% Part One : Analytical Solution

% Plot analytical solution
figure(1);
yline(0, '--');
grid on;

% Overdamped
plot(timeSeries, analyticalSolution(timeSeries, 2.5, 0.1, x_initial));
hold on;

% Critically Damped
plot(timeSeries, analyticalSolution(timeSeries, 2, 1, x_initial), 'g');
hold on;

% Underdamped
plot(timeSeries, analyticalSolution(timeSeries, 1, 5, x_initial), 'r');
title('Analytical Solution');
xlabel('Time');
ylabel('Position');
legend('Overdamped', 'Critically Damped', 'Underdamped');
hold off;

%% Part Two : Exploration of Analytical Solution
% Plot surface plot of gamma/k vs time to visualise how the 
% ratio (dampening) of the parameters affect the function

% Fix value of parameter k
kExplore = 1;

% Determine number of points required in discretization
noPoints = 100;

timeSeries = linspace(0, totalTime, noPoints);
gammaSeries = linspace(0, 2, noPoints);
ratioSeries = nan(size(gammaSeries));
positionSeries = nan(size(gammaSeries));

% build up ratioSeries
for i = 1:length(gammaSeries)
    % Calculate gamma using kExplore
    gamma = gammaSeries(i);

    % Calculate ratio using gamma and kExplore
    ratioSeries(i) = gamma.^2 ./ kExplore;

    % Calculate position based on timeSeries, gamma, kExplore and x_inital
    positionSeries(i, :) = analyticalSolution(timeSeries, gamma, kExplore, x_initial);
end

% Plot ratio of parameters vs position and time
figure(2);
surf(timeSeries, ratioSeries, positionSeries);

% Decorate surface plot
colorbar
title('Behaviour of anaytical solution')
xlabel('Gamma squared / k = 4');
ylabel('Time');
zlabel('Position');

% Adjust camera line of sight
view([-15 3 4]);

%% Part Three : Numberical Solution (USE EULERS OR ANOTHER METHOD)

% Convert 2nd order equation into a system of two coupled 1st order 
% equations using the variable suspension technique
%   Let xdot = v
%   Let vdot = xddot = -gamma * xdot - k * x
%
%   Since xdot = v, then...
%   d/dt [x; v] = [0 1; -k -gamma] * [x; v]

% TO REMOVE
k = 0.25; % frequency term (TO REMOVE)
gamma = 0.05; % dampening term (TO REMOVE)

% Negatively dampened (will converge at y = 0)
% Also known as an underdamped system
A = [0 1; -k -gamma];

% Use ode45 to numerically solve
[t, x] = ode45(@(t, x) A * x, timeSeries, x0);
position = x(:, 1);
velocity = x(:, 2);

%% Part Four : Plotting

% Plot Position vs Time
figure(3);
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

%% Part Five : Convergence

% Pick particular case, underdamped in this case
% Prove the numerical solution converges at the same value the analytical
% solution conveges to

% Analytical solution
analytical_position = analyticalSolution(timeSeries, 0.1, 3, x_initial);

% Numerical solution (grabbed from above but to be updated and turned into
% a function)
k = 0.25; % frequency term (TO REMOVE)
gamma = 0.05; % dampening term (TO REMOVE)

% Negatively dampened (will converge at y = 0)
% Also known as an underdamped system
A = [0 1; -k -gamma];

% Use ode45 to numerically solve
[t, x] = ode45(@(t, x) A * x, timeSeries, x0);
numerical_position = x(:, 1);

% Take normal of Absolute Error of analytical and numerical solution
% Reference: https://sutherland.che.utah.edu/wiki/index.php/Iteration_and_Convergence
normalAbsError = norm(abs(analytical_position - numerical_position));
fprintf('Normal of Absolute Error = %f\n', normalAbsError);

% Calculate Relative Error
normalRelError = 0; % TODO
fprintf('Normal of Absolute Error = %f\n', normalRelError);

% Expect Absolute and Relative error to indicate that the numerical
% solution is a valid model of the analytical solution
% TODO

%% Part Six : Analyis of error with varying step size

% Generate range of step sizes
stepSizes = linspace(0.001, 1);
errorAtStepSize = nan(length(stepSizes));

% Plot error vs step size
% We expect the error to be reduced as the step size is minimised
figure(4);
loglog(stepSizes, errorAtStepSize);
title('Error between Analytical Solution and Numerical Solutions vs Step Size');
xlabel('Step Size (log)');
ylabel('Normalised error');
legend('Error');

%% Part Seven : Fourier Analysis of numerical solution

% Do this for each case
% Compute a power spectrum graph of the numerical solution
% Plot (frequency, power_freqSpace);

% Determine the center frequency and FWHM
% Include these in the plot

%% Part Eight : Function Definitions

function position = analyticalSolution(timeSeries, gamma, k, x_init)  
    %   Derivation
    %   Let x = e^bt
    %   therefore... xdot = b * e^bt
    %   therefore... xddot = b^2 * e^bt
    %
    %   Plugging into the original PDE give us...
    %   -b^2 * e^bt - (gamma * b * e^bt) - ke^bt = 0
    %
    %   Pull e^bt out as common factor, this leaves us...
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
    
    % Define discriminant of characteristic equation
    discriminant = gamma.^2 - (4 .* k);
    
    % Define function in three cases based on the determinant of the roots
    % Reference: https://nrich.maths.org/11054
    if discriminant == 0 % critically damped
        % Solve for A and B constants
        A = x_init;
        B = x_init .* b_1;
        
        position = (A + B.*timeSeries) .* exp(b_1 .* timeSeries);
    elseif discriminant > 0 % overdamped
        % Solve for A and B constants
        A = (x_init .* b_2) ./ (b_2 - b_1);
        B = (x_init .* b_1) ./ (b_1 - b_2);
    
        position = (A .* exp(b_1 .* timeSeries)) + ...
            (B .* exp(b_2 .* timeSeries));

    else % underdamped if discriminant is less than 0
        % Separate real and imaginary parts of roots
        alpha = real(b_1);
        beta = imag(b_2);
        
        % Solve for A and B constants
        A = x_init;
        B = -x_init ./ beta;

        position = exp(alpha .* timeSeries) .* ...
            (A.*cos(beta.*timeSeries) + B.*sin(beta.*timeSeries));
    end
end
