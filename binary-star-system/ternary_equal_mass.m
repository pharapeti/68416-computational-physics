%% Set up workspace
clear; clc; close all

%% Set up Simulation
% Simulation Parameters
timestep = 60; % timestep (seconds)
totalTime = 60 * 60 * 24 * 15; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;
G = 6.67408 * 10.^-11;

%% Stability Check

x_series = linspace(-5e8, 5e8, 250);
y_series = linspace(-5e8, 5e8, 250);

% Determine whether system is stable in the following configurations
stableGrid = nan([length(x_series), length(y_series)]);

final_index = length(x_series) - 1;

for i = 1:length(x_series)
    for j = 1:length(y_series)

        % Create Simulation
        simulation = Simulation(timestep, totalTime, G);

        % Create two large body bodies
        star1_mass = 6 * 10.^24; % mass (kg)
        star1_initial_position = [0, 0]; % position (m)
        star1_initial_velocity = [0, -1e+3]; % velocity (m/s)
        star1_initial_acceleration = [0, 0]; % acceleration (m/s/s)
        simulation.createBody(star1_mass, star1_initial_position, ...
        star1_initial_velocity, star1_initial_acceleration);

        star2_mass = 6 * 10.^24; % mass (kg)
        star2_initial_position = [2e+8, 0]; % position (m)
        star2_initial_velocity = [0, 1e+3]; % velocity (m/s)
        star2_initial_acceleration = [0, 0]; % acceleration (m/s/s)
        simulation.createBody(star2_mass, star2_initial_position, ...
        star2_initial_velocity, star2_initial_acceleration);

        % Create Satellite
        x_position = x_series(:, i);
        y_position = y_series(:, j);
        
        satellite_mass = 6 * 10.^24; % mass (kg)
        satellite_initial_position = [x_position, y_position]; % position (m)
        satellite_initial_velocity = [-1e+3, 0]; % velocity (m/s)
        satellite_initial_acceleration = [0, 0]; % acceleration (m/s/s)
        simulation.createBody(satellite_mass, satellite_initial_position, ...
        satellite_initial_velocity, satellite_initial_acceleration);

        % Solve system over the duration of the simulation
        solveSystemNumerically(simulation);

        % Determine whether system is stable
        stableGrid(i, j) = 1;
        
        % Grab both stars and Satellite
        star1 = simulation.Bodies(1);
        star2 = simulation.Bodies(2);
        satellite = simulation.Bodies(3);

        % Calculate distances between bodies at end of simulation
        r1 = calculateDistanceBetweenBodies(final_index, star1, star2);
        r2 = calculateDistanceBetweenBodies(final_index, satellite, star1);
        r3 = calculateDistanceBetweenBodies(final_index, satellite, star2);
        
        % Check if all bodies haven't been ejected
        stable = 1;
        
        if (r1(1) > 1e+10 || r1(2) > 1e+10)
            stable = 0;
        end
        if (r2(1) > 1e+10 || r2(2) > 1e+10)
            stable = 0;
        end
        if (r3(1) > 1e+10 || r3(2) > 1e+10)
            stable = 0;
        end

        stableGrid(i, j) = stable;
    end
end


figure('NumberTitle', 'off', 'Name', 'Ejection map');
contourf(x_series, y_series, stableGrid);
title('Ejection map')
xlabel('X')
ylabel('Y')

%% Visualise System

% figure('NumberTitle', 'off', 'Name', 'Trajectory of Star 1, 2 and Satellite');
% plot(simulation.Bodies(1).Position.X, simulation.Bodies(1).Position.Y, 'r', ...
%      simulation.Bodies(2).Position.X, simulation.Bodies(2).Position.Y, 'b', ...
%      simulation.Bodies(3).Position.X, simulation.Bodies(3).Position.Y, 'g', ...
%      simulation.Barycenter.X, simulation.Barycenter.Y, '-')
% title('Trajectory of Star 1, 2 and Satellite');
% xlabel('Position X');
% ylabel('Position Y');
% legend('Star 1', 'Star 2', 'Satellite', 'Barycenter');
% 
% figure('NumberTitle', 'off', 'Name', 'Trajectory of Star 1');
% plot(simulation.Bodies(1).Position.X, simulation.Bodies(1).Position.Y);
% title('Tragectory of Star 1');
% xlabel('Position X');
% ylabel('Position Y');
% 
% figure('NumberTitle', 'off', 'Name', 'Trajectory of Star 2');
% plot(simulation.Bodies(2).Position.X, simulation.Bodies(2).Position.Y);
% title('Tragectory of Star 2');
% xlabel('Position X');
% ylabel('Position Y');
% 
% figure('NumberTitle', 'off', 'Name', 'Trajectory of Satellite');
% plot(simulation.Bodies(3).Position.X, simulation.Bodies(3).Position.Y);
% title('Tragectory of Satellite');
% xlabel('Position X');
% ylabel('Position Y');

%% Functions

% This function computes a numerical solution for a given dataset via the
% Verlet Integration's method.
function solveSystemNumerically(simulation)
    % Stars
    star1 = simulation.Bodies(1);
    star2 = simulation.Bodies(2);
    sat = simulation.Bodies(3);

    % Set up data structure to record barycenter across time
    barycenterHistory = Barycenter();
    barycenterHistory.X = nan(1, length(simulation.TimeSeries));
    barycenterHistory.Y = nan(1, length(simulation.TimeSeries));

    % Find initial barycenter
    initial_barycenter = barycenterFromOrigin(1, star1, star2);
    barycenterHistory.X(1) = initial_barycenter(1);
    barycenterHistory.Y(1) = initial_barycenter(2);

    % Inspired by http://www.astro.yale.edu/coppi/astro520/solving_differential_equation.pdf
    % Velocity Verlet Method
    % Numerically solve the position and velocity of the model system using
    % function arguments
    for i = 1:length(simulation.TimeSeries) - 1     

        % Calculate next position using previous position and second and
        % third order terms
        % Star 1
        star1.Position.X(i + 1) = star1.Position.X(i) + ...
            (simulation.TimeStep .* star1.Velocity.X(i)) + ...
            0.5 * star1.Acceleration.X(i) * simulation.TimeStep.^2;
        
        star1.Position.Y(i + 1) = star1.Position.Y(i) + ...
            (simulation.TimeStep .* star1.Velocity.Y(i)) + ...
            0.5 * star1.Acceleration.Y(i) * simulation.TimeStep.^2;
        
        % Star 2
        star2.Position.X(i + 1) = star2.Position.X(i) + ...
            (simulation.TimeStep .* star2.Velocity.X(i)) + ...
            0.5 * star2.Acceleration.X(i) * simulation.TimeStep.^2;
        
        star2.Position.Y(i + 1) = star2.Position.Y(i) + ...
            (simulation.TimeStep .* star2.Velocity.Y(i)) + ...
            0.5 * star2.Acceleration.Y(i) * simulation.TimeStep.^2;
        
        % Satellite
        sat.Position.X(i + 1) = sat.Position.X(i) + ...
            (simulation.TimeStep .* sat.Velocity.X(i)) + ...
            0.5 * sat.Acceleration.X(i) * simulation.TimeStep.^2;
        
        sat.Position.Y(i + 1) = sat.Position.Y(i) + ...
            (simulation.TimeStep .* sat.Velocity.Y(i)) + ...
            0.5 * sat.Acceleration.Y(i) * simulation.TimeStep.^2;
        
        % Forces between Earth and Moon
        d_earth_moon = calculateDistanceBetweenBodies(i + 1, star1, star2);
        d_length = norm(d_earth_moon);
        unit_d = d_earth_moon ./ d_length;
        
        force_earth_moon = (simulation.G * star1.Mass * star2.Mass ./ d_length .^ 2) .* unit_d;

        % Forces between Earth and Satellite
        d_sat_earth = calculateDistanceBetweenBodies(i + 1, star1, sat);
        d_length = norm(d_sat_earth);
        unit_d = d_sat_earth ./ d_length;
        
        force_sat_earth = (simulation.G * star1.Mass * sat.Mass ./ d_length .^ 2) .* unit_d;
        
        % Forces between Moon and Satellite
        d_sat_moon = calculateDistanceBetweenBodies(i + 1, star2, sat);
        d_length = norm(d_sat_moon);
        unit_d = d_sat_moon ./ d_length;

        force_sat_moon = (simulation.G * star2.Mass * sat.Mass ./ d_length .^ 2) .* unit_d;

        % Using Newton's Third Law, calculate equal and opposite forces
        force_moon_earth = -force_earth_moon;
        force_earth_sat = -force_sat_earth;
        force_moon_sat = -force_sat_moon;

        % Calculate sum of forces of each body
        sum_force_earth = force_earth_moon + force_earth_sat;
        sum_force_moon = force_moon_earth + force_moon_sat;
        sum_force_sat = -(force_sat_earth + force_sat_moon);

        % Using Newton's Second Law on each mass, determine accelerations 
        % on each star in (x, y) directions due to forces
        star1.Acceleration.X(i + 1) = sum_force_earth(1) ./ star1.Mass;
        star1.Acceleration.Y(i + 1) = sum_force_earth(2) ./ star1.Mass;
        
        star2.Acceleration.X(i + 1) = sum_force_moon(1) ./ star2.Mass;
        star2.Acceleration.Y(i + 1) = sum_force_moon(2) ./ star2.Mass;
        
        sat.Acceleration.X(i + 1) = sum_force_sat(1) ./ sat.Mass;
        sat.Acceleration.Y(i + 1) = sum_force_sat(2) ./ sat.Mass;
        
        % Calculate new velocity based on previous velocity and third order
        % term
        % Star 1
        star1.Velocity.X(i + 1) = star1.Velocity.X(i) + ...
            0.5 * (star1.Acceleration.X(i) + ...
            star1.Acceleration.X(i + 1)) * simulation.TimeStep;

        star1.Velocity.Y(i + 1) = star1.Velocity.Y(i) + ...
            0.5 * (star1.Acceleration.Y(i) + ...
            star1.Acceleration.Y(i + 1)) * simulation.TimeStep;
        
        % Star 2
        star2.Velocity.X(i + 1) = star2.Velocity.X(i) + ...
            0.5 * (star2.Acceleration.X(i) + ...
            star2.Acceleration.X(i + 1)) * simulation.TimeStep;
        
        star2.Velocity.Y(i + 1) = star2.Velocity.Y(i) + ...
            0.5 * (star2.Acceleration.Y(i) + ...
            star2.Acceleration.Y(i + 1)) * simulation.TimeStep;
        
        % Satellite
        sat.Velocity.X(i + 1) = sat.Velocity.X(i) + ...
            0.5 * (sat.Acceleration.X(i) + ...
            sat.Acceleration.X(i + 1)) * simulation.TimeStep;
        
        sat.Velocity.Y(i + 1) = sat.Velocity.Y(i) + ...
            0.5 * (sat.Acceleration.Y(i) + ...
            sat.Acceleration.Y(i + 1)) * simulation.TimeStep;

        % Calculate Barycenter of system
        barycenter = barycenterFromOrigin(i, star1, star2);
        barycenterHistory.X(i) = barycenter(1);
        barycenterHistory.Y(i) = barycenter(2);
    end

    % Persist changes to simulation - due to inability to edit by reference
    simulation.Bodies(1) = star1;
    simulation.Bodies(2) = star2;
    simulation.Bodies(3) = sat;
    simulation.Barycenter.X = barycenterHistory.X;
    simulation.Barycenter.Y = barycenterHistory.Y;
end

% Returns 1x2 array of distances in (x, y) respectively
function r = calculateDistanceBetweenBodies(i, bodyA, bodyB)
    positionA = [bodyA.Position.X(i), bodyA.Position.Y(i)];
    positionB = [bodyB.Position.X(i), bodyB.Position.Y(i)];

    r = positionB - positionA;
end

function origin_to_barycenter = barycenterFromOrigin(i, bodyA, bodyB)
    % From perspective of star 1
    distance_between_bodies = calculateDistanceBetweenBodies(i, bodyA, bodyB);
    mass_ratio = bodyA.Mass ./ bodyB.Mass;

    % We should calculate the barycenter in respect to the (x, y) 
    % coordinate origin to allow our measurements to be agnostic to a
    % specific body in the system.

    % We can acheive this by imagine the vector Rb to be the vector from
    % the coordinate origin (x, y = [0, 0]) to the position of a particular
    % body.
    origin_to_star_A = [bodyA.Position.X(i), bodyA.Position.Y(i)];

    % Then we can imagine another vector Rc to be the vector from the
    % position of a chosen body in the system to the position of the
    % system's barycenter.
    star_A_to_barycenter = distance_between_bodies ./ (1 + mass_ratio);
    
    % Using these two vectors (Rb and Rc), we can compose a resultant
    % vector Rbc which is the sum of Rb and Rc. The vector addition will
    % result in Rbc being to vector from the coordinate origin (xy = [0,
    % 0]) to the position of the barycenter of the system.
    origin_to_barycenter = origin_to_star_A + star_A_to_barycenter;
end
