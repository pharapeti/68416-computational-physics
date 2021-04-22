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
star1_initial_position = [100, 100];
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

figure(5);
plot(simulation.Barycenter(:, 1), simulation.Barycenter(:, 2));
grid on;
title('Barycenter (x, y)');

%% Functions

% This function computes a numerical solution for a given dataset via the
% Euler's method.
function solveSystemNumerically(simulation)
    % Stars
    star1 = simulation.Bodies(1);
    star2 = simulation.Bodies(2);
    
    % Set up data structure to record barycenter across time
    barycenterHistory = nan(length(simulation.TimeSeries), 2);

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

        barycenterHistory(i, :) = barycenterFromOrigin(i, star1, star2);
    end

    % Persist changes to simulation - due to inability to edit by reference
    simulation.Bodies(1) = star1;
    simulation.Bodies(2) = star2;
    simulation.Barycenter = barycenterHistory;
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

function origin_to_barycenter = barycenterFromOrigin(i, star1, star2)
    % From perspective of star 1
    distance_between_stars = calculateDistanceBetweenStars(i, star1, star2);
    mass_ratio = star1.Mass ./ star2.Mass;

    % We should calculate the barycenter in respect to the (x, y) 
    % coordinate origin to allow our measurements to be agnostic to a
    % specific body in the system.
    
    % We can acheive this by imagine the vector Rb to be the vector from
    % the coordinate origin (x, y = [0, 0]) to the position of a particular
    % body.
    origin_to_star_1 = star1.Position(i);
    
    % Then we can imagine another vector Rc to be the vector from the
    % position of a chosen body in the system to the position of the
    % system's barycenter.
    star_1_to_barycenter = distance_between_stars ./ (1 + mass_ratio);
    
    % Using these two vectors (Rb and Rc), we can compose a resultant
    % vector Rbc which is the sum of Rb and Rc. The vector addition will
    % result in Rbc being to vector from the coordinate origin (xy = [0,
    % 0]) to the position of the barycenter of the system.
    origin_to_barycenter = origin_to_star_1 + star_1_to_barycenter;
end
