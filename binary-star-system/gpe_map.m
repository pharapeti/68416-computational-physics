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
star1_initial_acceleration = [0, 0]; % acceleration (m/s/s)
simulation.createBody(star1_mass, star1_initial_position, star1_initial_velocity, star1_initial_acceleration);

% Create Star 2
star2_mass = 7.34767309 * 10.^22; % mass (kg)
star2_initial_position = [384400 * 10.^3, 0]; % position (m)
star2_initial_velocity = [0, 1.022 * 10.^3]; % velocity (m/s)
star2_initial_acceleration = [0, 0]; % acceleration (m/s/s)
simulation.createBody(star2_mass, star2_initial_position, star2_initial_velocity, star2_initial_acceleration);

% Set up Star 2 with velocities, such that it orbits Star 1
% simulation.Bodies(2).Velocity(1:0) = initialVelocityRequiredForOrbit(simulation, simulation.Bodies(1), simulation.Bodies(1));

%% Solve System
solveSystemNumerically(simulation);

%% Generate contour plot of Effective Graviational Potential Energy

% Step through values of x and y for a instance of  a simulation and 
% calculate the Effective Graviational Potential Energy at each point

% For example, take the first iteration of the Earth-Moon system
x_series = linspace(-5e8, 5e8, 1000);
y_series = linspace(-5e8, 5e8, 1000);

% Calculate Effective GPE over the x, y domain
gpe = nan([length(x_series), length(y_series)]);
for i = 1:length(x_series)
    for j = 1:length(y_series)
        x_position = x_series(:, i);
        y_position = y_series(:, j);

        % Calculate effective GPE at this point
        gpe(i, j) = calculateEffectiveGPE(simulation, x_position, y_position);
    end
end

% Generate logarithmic contour plot of Effective GPE
figure(15);
log_gpe = log10(abs(gpe));
contourf(x_series, y_series, log_gpe.');
colorbar;
grid on;
title('Effective Graviational Potential Energy');
xlabel('Y');
ylabel('X');
hold on;

% Mark position of Earth
x = simulation.Bodies(1).Position.X(1);
y = simulation.Bodies(1).Position.Y(1);
plot(x, y, '+g', 'LineWidth', 4, 'MarkerSize', 12);

% Mark position of Moon
x = simulation.Bodies(2).Position.X(1);
y = simulation.Bodies(2).Position.Y(1);
plot(x, y, '*', 'LineWidth', 4, 'MarkerSize', 12);

% Find Lagrange Points (maximums)
maximum = max(max(gpe));
[x_index, y_index] = find(gpe == maximum);
text(x_series(x_index(2)), y_series(y_index(2)), 'L4', 'Color', 'white', 'FontSize', 16)
text(x_series(x_index(1)), y_series(y_index(1)), 'L5', 'Color', 'white', 'FontSize', 16)

% Find Lagrange Points (minimums)
minimum = min(min(gpe));
[x_index, y_index] = find(gpe == minimum);
text(x_series(x_index(2)), y_series(y_index(2)), 'L1', 'Color', 'black', 'FontSize', 16)
text(x_series(x_index(1)), y_series(y_index(1)), 'L2', 'Color', 'black', 'FontSize', 16)

hold off;

%% Functions

% This function computes a numerical solution for a given dataset via the
% Verlet Integration's method.
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

        % Calculate new acceleration based on new position
        % Calculate distance between star 1 and 2
        distance = calculateDistanceBetweenBodies(i + 1, star1, star2);
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
        star1.Acceleration.X(i + 1) = forceXOn1 ./ star1.Mass;
        star1.Acceleration.Y(i + 1) = forceYOn1 ./ star1.Mass;
        
        star2.Acceleration.X(i + 1) = forceXOn2 ./ star2.Mass;
        star2.Acceleration.Y(i + 1) = forceYOn2 ./ star2.Mass;
        
        % Calculate new velocity based on previous velocity and third order
        % term
        % Star 1
        star1.Velocity.X(i + 1) = star1.Velocity.X(i) + ...
            0.5 * (star1.Acceleration.X(i) + star1.Acceleration.X(i + 1)) * simulation.TimeStep;

        star1.Velocity.Y(i + 1) = star1.Velocity.Y(i) + ...
            0.5 * (star1.Acceleration.Y(i) + star1.Acceleration.Y(i + 1)) * simulation.TimeStep;
        
        % Star 2
        star2.Velocity.X(i + 1) = star2.Velocity.X(i) + ...
            0.5 * (star2.Acceleration.X(i) + star2.Acceleration.X(i + 1)) * simulation.TimeStep;
        
        star2.Velocity.Y(i + 1) = star2.Velocity.Y(i) + ...
            0.5 * (star2.Acceleration.Y(i) + star2.Acceleration.Y(i + 1)) * simulation.TimeStep;

        % Calculate Barycenter of system
        barycenter = barycenterFromOrigin(i, star1, star2);
        barycenterHistory.X(i) = barycenter(1);
        barycenterHistory.Y(i) = barycenter(2);
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

function GPE = calculateEffectiveGPE(sim, x_position, y_position)
    current_position = [x_position, y_position];
    total_mass_to_distance_ratio = 0;

    for i = 1:length(sim.Bodies) 
        body = sim.Bodies(i);
        body_position = [body.Position.X(1), body.Position.Y(1)];

        % Find distance between current position and body
        distance_vector = body_position - current_position;
        distance_vector_length = norm(distance_vector);

        mass_to_distance_ratio = body.Mass ./ distance_vector_length;
        total_mass_to_distance_ratio = total_mass_to_distance_ratio + ...
            mass_to_distance_ratio;
    end
    
    % Ref http://www.physics.usyd.edu.au/~helenj/SeniorAstro/lecture11.pdf
    % Divide both sides by negligible m to give you effective potential
    % energy per unit mass
    vector_between_bodies = calculateDistanceBetweenBodies(1, sim.Bodies(1), sim.Bodies(2));
    distance_between_bodies = norm(vector_between_bodies);

    total_mass = sim.Bodies(1).Mass + sim.Bodies(2).Mass;
    omega_squared = (sim.G * total_mass) ./ distance_between_bodies.^3;

    barycenter_location = barycenterFromOrigin(1, sim.Bodies(1), sim.Bodies(2));
    r = norm(barycenter_location - current_position);
    centripedal_potential = 0.5 .* omega_squared * r^2;

    GPE = (-sim.G .* total_mass_to_distance_ratio) - centripedal_potential;
end
