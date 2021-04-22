%% Set up workspace
clear; clc; close all

%% Set up Simulation
% Simulation Parameters
timestep = 1; % timestep (seconds)
totalTime = 10e2; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;
G = 6.67408 * 10.^-11;

% Create Simulation
simulation = Simulation(timestep, totalTime, G);

% Create Star 1
star1_mass = 15;
star1_initial_position = [1, 1];
star1_initial_velocity = [0, 0];
simulation.createBody(star1_mass, star1_initial_position, star1_initial_velocity);

% Create Star 2
star2_mass = 10;
star2_initial_position = [0, 0];
star2_initial_velocity = [0, 0];
simulation.createBody(star2_mass, star2_initial_position, star2_initial_velocity);

% Set up Star 2 with velocities, such that it orbits Star 1
% simulation.Bodies(2).Velocity(1:0) = initialVelocityRequiredForOrbit(simulation, simulation.Bodies(1), simulation.Bodies(1));

%% Solve System
solveSystemNumerically(simulation);

%% Visualise System
 
figure(1);
plot3(simulation.TimeSeries, simulation.Bodies(1).Position(:, 1), simulation.Bodies(1).Position(:, 2));
grid on;
title('Star 1 - Position X,Y vs Time');

figure(2);
plot3(simulation.TimeSeries, simulation.Bodies(2).Position(:, 1), simulation.Bodies(2).Position(:, 2));
grid on;
title('Star 2 - Position X,Y vs Time');

figure(3);
plot(simulation.Bodies(1).Position(:, 1), simulation.Bodies(1).Position(:, 2));
grid on;
title('Star 1 - Position X vs Position Y');

figure(4);
plot(simulation.Bodies(2).Position(:, 1), simulation.Bodies(2).Position(:, 2));
grid on;
title('Star 2 - Position X vs Position Y');

%% Functions

% This function computes a numerical solution for a given dataset via the
% Euler's method.
function solveSystemNumerically(simulation)
    % Stars
    star1 = simulation.Bodies(1);
    star2 = simulation.Bodies(2);

    % Gravitational Factor
    g_factor = simulation.G * star1.Mass * star2.Mass;

    % Inspired by http://www.astro.yale.edu/coppi/astro520/solving_differential_equation.pdf
    % Numerically solve the position and velocity of the model system using
    % function arguments
    for i = 1:length(simulation.TimeSeries) - 1
        % Calculate distance between star 1 and 2
        distance = calculateDistanceBetweenStars(i, star1, star2);
        distance_vector_length = norm(distance);
        unit_distance = distance ./ distance_vector_length;

        % Using distance, calculate force on each star
        forceOn1 = (g_factor ./ (distance_vector_length .^ 2)) .* unit_distance;

        % Split forces into (x, y) directions
        forceXOn1 = forceOn1(1);
        forceYOn1 = forceOn1(2);

        % Newton's Second Law - equal and opposite reaction force
        forceXOn2 = - forceXOn1;
        forceYOn2 = - forceYOn1;

        % Using force on each mass, determine accelerations on each star in
        % (x, y) directions
        accelerationXOf1 = forceXOn1 ./ star1.Mass;
        accelerationYOf1 = forceYOn1 ./ star1.Mass;

        accelerationXOf2 = forceXOn2 ./ star2.Mass;
        accelerationYOf2 = forceYOn2 ./ star2.Mass;
        
        % Update Velocities
        % Star 1
        star1.Velocity(i + 1, 1) = star1.Velocity(i, 1) + (simulation.TimeStep .* accelerationXOf1);
        star1.Velocity(i + 1, 2) = star1.Velocity(i, 2) + (simulation.TimeStep .* accelerationYOf1);
        
        % Star 2
        star2.Velocity(i + 1, 1) = star2.Velocity(i, 1) + (simulation.TimeStep .* accelerationXOf2);
        star2.Velocity(i + 1, 2) = star2.Velocity(i, 2) + (simulation.TimeStep .* accelerationYOf2);

        % Update Positions
        % Star 1
        star1.Position(i + 1, 1) = star1.Position(i, 1) + (simulation.TimeStep .* star1.Velocity(i, 1));
        star1.Position(i + 1, 2) = star1.Position(i, 2) + (simulation.TimeStep .* star1.Velocity(i, 2));
        
        % Star 2
        star2.Position(i + 1, 1) = star2.Position(i, 1) + (simulation.TimeStep .* star2.Velocity(i, 1));
        star2.Position(i + 1, 2) = star2.Position(i, 2) + (simulation.TimeStep .* star2.Velocity(i, 2));
    end
    
    % Persist changes to simulation - due to inability to edit by reference
    simulation.Bodies(1) = star1;
    simulation.Bodies(2) = star2;
end

% Returns 1x2 array of distances in (x, y) respectively
function r = calculateDistanceBetweenStars(i, star_A, star_B)
    starA = star_A.Position(i, :);
    starB = star_B.Position(i, :);

    r = starB - starA;
end

% Returns the orbital velocity required for star B to orbit star A
% Only holds true when the mass of Star A is much greater than that of Star
% B
function velocity = initialVelocityRequiredForOrbit(simulation, star_A, star_B)
    r = abs(calculateDistanceBetweenStars(1, star_A, star_B));
    velocity = sqrt((simulation.G * star_A.Mass) ./ r);
end
