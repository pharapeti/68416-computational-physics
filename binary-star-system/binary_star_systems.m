%% Part Zero : Setup

% Clear existing workspace
clear; clc; close all

% Setup parameters
timestep = 0.01; % timestep (seconds)
totalTime = 10; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;

% Constants
G = 6.67 * 10.^-11;

% Star 1
star1_mass = 10;
star1_initial_position = [5, 5];
star1_initial_velocity = [0, 0];
star1 = Body(star1_mass, star1_initial_position, star1_initial_velocity);

% Star 2
star2_mass = 20;
star2_initial_position = [0, 0];
star2_initial_velocity = [0, 0];
star2 = Body(star2_mass, star2_initial_position, star2_initial_velocity);

starsInSystem = [star1, star2];

solution = generateNumericalSolution(starsInSystem, timeSeries, timestep);
star_1_motion = solution(:, 1:4);
star_2_motion = solution(:, 5: 8);

figure(1);
plot3(timeSeries, star_1_motion(:, 1), star_1_motion(:, 2));
grid on;
title('Star 1 - X,Y vs Time');

figure(2);
plot3(timeSeries, star_2_motion(:, 1), star_2_motion(:, 2));
grid on;
title('Star 2 - X,Y vs Time');
% 
% subplot(1, 8, 2);
% plot(timeSeries, positionXOf2);
% title('PosX of 2');
% 
% subplot(1, 8, 3);
% plot(timeSeries, positionYOf1);
% title('PosY of 1');
% 
% subplot(1, 8, 4);
% plot(timeSeries, positionYOf2);
% title('PosY of 2');
% 
% subplot(1, 8, 5);
% plot(timeSeries, velocityXOf1);
% title('VelX of 1');
% 
% subplot(1, 8, 6);
% plot(timeSeries, velocityXOf2);
% title('VelX of 2');
% 
% subplot(1, 8, 7);
% plot(timeSeries, velocityYOf1);
% title('VelY of 1');
% 
% subplot(1, 8, 8);
% plot(timeSeries, velocityYOf2);
% title('VelY of 2');

% This function computes a numerical solution for a given dataset via the
% Euler's method.
function starMotionArray = generateNumericalSolution(starsInSystem, timeSeries, stepSize)
    % Stars
    star1 = starsInSystem(1);
    star2 = starsInSystem(2);

    % Generate data structure for each star
    % Rows represent time steps
    star1_Motion = nan(length(timeSeries), 4);
    star2_Motion = nan(length(timeSeries), 4);
    
    % Column 1 is position in x direction
    star1_Motion(1, 1) = star1.Position(1);
    star2_Motion(1, 1) = star2.Position(1);
    
    % Column 2 is position in y direction
    star1_Motion(1, 2) = star1.Position(2);
    star2_Motion(1, 2) = star2.Position(2);
    
    % Column 3 is velocity in x direction
    star1_Motion(1, 3) = star1.Velocity(1);
    star2_Motion(1, 3) = star2.Velocity(1);
    
    % Column 4 is velocity in y direction
    star1_Motion(1, 4) = star1.Velocity(2);
    star2_Motion(1, 4) = star2.Velocity(2);

    % Stars
    star1 = starsInSystem(1);
    star2 = starsInSystem(2);
    
    % Gravitational Factor
    G = 6.67408 * 10.^-11;
    g_factor = G * star1.Mass * star2.Mass;

    % Inspired by http://www.astro.yale.edu/coppi/astro520/solving_differential_equation.pdf
    % Numerically solve the position and velocity of the model system using
    % function arguments
    for i = 1:length(timeSeries) - 1
        % Calculate distance between star 1 and 2
        vector_r = calculateDistanceBetweenStars(star1, star2);
        vector_r_length = norm(vector_r);
        unit_vector_r = vector_r ./ vector_r_length;

        % Using distance, calculate force on each star
        forceOn1 = (g_factor ./ (vector_r_length .^ 2)) .* unit_vector_r;
        
        % Split forces into (x, y) directions
        forceXOn1 = forceOn1(1);
        forceYOn1 = forceOn1(2);

        % Newton's Second Law
        % Equal and opposite reaction force
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
        star1_Motion(i + 1, 3) = star1_Motion(i, 3) + (stepSize .* accelerationXOf1);
        star1_Motion(i + 1, 4) = star1_Motion(i, 4) + (stepSize .* accelerationYOf1);
        
        % Star 2
        star2_Motion(i + 1, 3) = star2_Motion(i, 3) + (stepSize .* accelerationXOf2);
        star2_Motion(i + 1, 4) = star2_Motion(i, 4) + (stepSize .* accelerationYOf2);

        % Update Positions
        % Star 1
        star1_Motion(i + 1, 1) = star1_Motion(i, 1) + (stepSize .* star1_Motion(i, 3));
        star1_Motion(i + 1, 2) = star1_Motion(i, 2) + (stepSize .* star1_Motion(i, 4));
        
        % Star 2
        star2_Motion(i + 1, 1) = star2_Motion(i, 1) + (stepSize .* star2_Motion(i, 3));
        star2_Motion(i + 1, 2) = star2_Motion(i, 2) + (stepSize .* star2_Motion(i, 4));
    end

    % Package and return motions of stars
    starMotionArray = [star1_Motion, star2_Motion];
end

% Returns 1x2 array of distances in (x, y) respectively
function r = calculateDistanceBetweenStars(star_A, star_B)
    % TODO just reference these directly instead of assigned to temporary
    % variable
    starA_x = star_A.Position(1);
    starA_y = star_A.Position(2);
    
    starB_x = star_B.Position(1);
    starB_y = star_B.Position(2);

    r = [starB_x - starA_x, starB_y - starA_y];
end

