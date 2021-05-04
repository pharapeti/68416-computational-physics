%% Set up workspace
clear; clc; close all

%% Set up Simulation
% Simulation Parameters
timestep = 60; % timestep (seconds)
totalTime = 60 * 60 * 24 * 27; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;
G = 6.67408 * 10.^-11;

% Create Simulation
simulation = Simulation(timestep, totalTime, G);

% Create Star 1
star1_mass = 5.97219 * 10.^24; % mass (kg)

% Earth's COM -> Barycenter = 4675km
star1_initial_position = [-4675 * 10.^3, 0]; % position (m)

% Earth orbit's the Earth-Moon Barycenter at 45km/h = 12.5m/ss
star1_initial_velocity = [0, 12.5]; % velocity (m/s)
simulation.createBody(star1_mass, star1_initial_position, star1_initial_velocity);

% Create Star 2
star2_mass = 7.34767309 * 10.^22; % mass (kg)
star2_initial_position = [384400 * 10.^3, 0]; % position (m)
star2_initial_velocity = [0, 1.022 * 10.^3]; % velocity (m/s)
simulation.createBody(star2_mass, star2_initial_position, star2_initial_velocity);

% Set up Star 2 with velocities, such that it orbits Star 1
% simulation.Bodies(2).Velocity(1:0) = initialVelocityRequiredForOrbit(simulation, simulation.Bodies(1), simulation.Bodies(1));

%% Solve System
solveSystemNumerically(simulation);

%% Visualise System

figure(1);
plot(simulation.Bodies(1).Position.X, simulation.Bodies(1).Position.Y, 'r', ...
     simulation.Bodies(2).Position.X, simulation.Bodies(2).Position.Y, 'b', ...
     simulation.Barycenter.X, simulation.Barycenter.Y, '-')
title('Orbit of Moon about the Earth');
xlabel('Position X');
ylabel('Position Y');
legend('Earth', 'Moon', 'Barycenter');

figure(2);
plot(simulation.Bodies(1).Position.X, simulation.Bodies(1).Position.Y);
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('Earth');
xlabel('Position X');
ylabel('Position Y');

figure(3);
plot(simulation.Bodies(2).Position.X, simulation.Bodies(2).Position.Y);
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('Moon');
xlabel('Position X');
ylabel('Position Y');

figure(4);
plot(simulation.TimeSeries, simulation.Bodies(1).PE);
hold on;
plot(simulation.TimeSeries, simulation.Bodies(1).KE, 'r');
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('KE and PE of Earth');
legend('PE', 'KE');
xlabel('Time (s)');
ylabel('Energy (j)');

figure(5);
plot(simulation.TimeSeries, simulation.Bodies(2).PE);
hold on;
plot(simulation.TimeSeries, simulation.Bodies(2).KE, 'r');
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('KE and PE of Moon');
legend('PE', 'KE');
xlabel('Time (s)');
ylabel('Energy (j)');

figure(6);
plot(simulation.TimeSeries, simulation.Bodies(1).PE);
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('PE of Moon');
legend('PE');
xlabel('Time (s)');
ylabel('Energy (j)');

figure(7);
plot(simulation.TimeSeries, simulation.Bodies(2).PE);
% xlim([-10.^6, 10.^6]);
% ylim([-10.^6, 10.^6]);
title('PE of Earth');
legend('PE');
xlabel('Time (s)');
ylabel('Energy (j)');

%% Validate Results with Energy Method

% Determine last index of simulation timeseries
final_time_index = length(simulation.TimeSeries) - 1;

% Compare total energy of system at first timestep to the total energy of
% the system at the final timestep
initial_energy = totalEnergyOfSimulationAtI(1, simulation);
final_energy = totalEnergyOfSimulationAtI(final_time_index, simulation);
delta_E = abs(final_energy - initial_energy);

% Find the relative error % between the two
relative_error = delta_E ./ final_energy;

% Do this for a range of step sizes to prove that minimising the step size
% also minimises the relative error

%% Validate Results with final position of the Moon convergence

step_sizes = logspace(2, 4);
final_positions_of_moon = nan(length(step_sizes), 2);

for i = 1:length(step_sizes)
    % Setup simulation at a given timestep
    sim = Simulation(step_sizes(i), totalTime, G);
    sim.createBody(star1_mass, star1_initial_position, star1_initial_velocity);
    sim.createBody(star2_mass, star2_initial_position, star2_initial_velocity);

    % Solve system numerically
    solveSystemNumerically(sim);

    % Determine final position of the Moon
    final_pos_x = simulation.Bodies(2).Position.X(final_time_index);
    final_pos_y = simulation.Bodies(2).Position.Y(final_time_index);
    final_positions_of_moon(i, :) = [final_pos_x, final_pos_y];
end

% Plot Step size vs Final Position of Moon
figure(8)
loglog(step_sizes, final_positions_of_moon(:, 1));
title('Position X of Moon at end of simulation vs Step Size');
xlabel('Step Size (log)');
ylabel('Position of Moon X');
legend('Position of Moon X');

% Plot Step size vs Final Position of Moon
figure(9)
loglog(step_sizes, final_positions_of_moon(:, 2));
title('Position Y of Moon at end of simulation vs Step Size');
xlabel('Step Size (log)');
ylabel('Position of Moon Y');
legend('Position of Moon Y');

% Do this for a range of step sizes to prove that minimising the step size
% also minimises the relative error

%% Functions

% This function computes a numerical solution for a given dataset via the
% Euler's method.
function solveSystemNumerically(simulation)
    % Stars
    star1 = simulation.Bodies(1);
    star2 = simulation.Bodies(2);

    % Set up data structure to record barycenter across time
    barycenterHistory = Barycenter();
    barycenterHistory.X = nan(1, length(simulation.TimeSeries));
    barycenterHistory.Y = nan(1, length(simulation.TimeSeries));

    % Gravitational Factor
    g_factor = simulation.G * star1.Mass * star2.Mass;

    % Find initial barycenter
    iniital_barycenter = barycenterFromOrigin(1, star1, star2);
    barycenterHistory.X(1) = iniital_barycenter(1);
    barycenterHistory.Y(1) = iniital_barycenter(2);

    % Inspired by http://www.astro.yale.edu/coppi/astro520/solving_differential_equation.pdf
    % Numerically solve the position and velocity of the model system using
    % function arguments
    for i = 1:length(simulation.TimeSeries) - 1
        % Calculate distance between star 1 and 2
        distance = calculateDistanceBetweenBodies(i, star1, star2);
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
        star1.Velocity.X(i + 1) = star1.Velocity.X(i) + (simulation.TimeStep .* accelerationXOf1);
        star1.Velocity.Y(i + 1) = star1.Velocity.Y(i) + (simulation.TimeStep .* accelerationYOf1);
        
        % Star 2
        star2.Velocity.X(i + 1) = star2.Velocity.X(i) + (simulation.TimeStep .* accelerationXOf2);
        star2.Velocity.Y(i + 1) = star2.Velocity.Y(i) + (simulation.TimeStep .* accelerationYOf2);

        % Update Positions
        % Star 1
        star1.Position.X(i + 1) = star1.Position.X(i) + (simulation.TimeStep .* star1.Velocity.X(i));
        star1.Position.Y(i + 1) = star1.Position.Y(i) + (simulation.TimeStep .* star1.Velocity.Y(i));
        
        % Star 2
        star2.Position.X(i + 1) = star2.Position.X(i) + (simulation.TimeStep .* star2.Velocity.X(i));
        star2.Position.Y(i + 1) = star2.Position.Y(i) + (simulation.TimeStep .* star2.Velocity.Y(i));

        barycenter = barycenterFromOrigin(i, star1, star2);
        barycenterHistory.X(i) = barycenter(1);
        barycenterHistory.Y(i) = barycenter(2);

        % Calculate Kinetic Energy at this timestep
        star1.KE(i) = calculateKE(i, star1);
        star2.KE(i) = calculateKE(i, star2);

        % Calculate Graviational Potential Energy at this timestep
        star1.PE(i) = calculatePE(i, star1, star2, simulation.G);
        star2.PE(i) = calculatePE(i, star2, star1, simulation.G);
    end

    % Persist changes to simulation - due to inability to edit by reference
    simulation.Bodies(1) = star1;
    simulation.Bodies(2) = star2;
    simulation.Barycenter.X = barycenterHistory.X;
    simulation.Barycenter.Y = barycenterHistory.Y;
end

% Returns 1x2 array of distances in (x, y) respectively
function r = calculateDistanceBetweenBodies(i, bodyA, bodyB)
    starA = [bodyA.Position.X(i), bodyA.Position.Y(i)];
    starB = [bodyB.Position.X(i), bodyB.Position.Y(i)];

    r = starB - starA;
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

% Calculates the Kinetic Energy for a given star
function KE = calculateKE(i, body)
    KE = 0.5 * body.Mass * ...
        (body.Position.X(i) .^2 + body.Position.Y(i) .^2);
end

% Calculates the Gravitational Potential Energy between two stars
function PE = calculatePE(i, bodyA, bodyB, G)
    r = calculateDistanceBetweenBodies(i, bodyA, bodyB);
    r_length = norm(r);
    PE = (-G .* bodyA.Mass .* bodyB.Mass) ./ r_length;
end

function E = totalEnergyOfSimulationAtI(i, simulation)
    E = 0;

    for j = 1:length(simulation.Bodies) - 1
        PE = simulation.Bodies(j).PE(i);
        KE = simulation.Bodies(j).KE(i);
        E = E + PE + KE;
    end
end
