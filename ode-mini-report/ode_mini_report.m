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
% ratio (dampening) of the parameters affect the function

%% Part Four : Convergence

% Pick particular case
% Prove the numerical solution converges at the same value the analytical
% solution conveges to

% Analytical solution
analytical_position = analyticalSolution(timeSeries, 0.1, 3, x_initial);

% Solve for constants
c = analytical_position.' \ position;

% Sum squared error
S = sum((analytical_position - position) .^ 2);

% Mean squared error
% ... m = number of data points = length(timeSeries)
% ... n = number of coefficients = length(c)
% ... therefore m - n
m = length(timeSeries);
n = length(c);
meanSquaredError = sqrt(S / (m - n));
fprintf('Mean Squared Error = %f\n', meanSquaredError);

% Generate Covariance Matrix from diagonal of inverse conjugate matrix of A
%covarianceMatrix = diag(inv(analytical_position.' * analytical_position));

% Uncertainty Matrix
%uncertainty = meanSquaredError .* sqrt(covarianceMatrix);

%% Part Five : Error between analytical and numerical solution

% Pick particular case
% Produce a plot of the error between the analytical and numerical
% solution as as the step size of increased/decreased
% Plot should be plot(stepsize, error)

%% Part Six : Fourier Analysis of numerical solution

% Do this for each case
% Compute a power spectrum graph of the numerical solution
% Plot (frequency, power_freqSpace);

% Determine the center frequency and FWHM
% Include these in the plot

%% Part Seven : Function Definitions

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
