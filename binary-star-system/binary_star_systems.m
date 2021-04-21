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

% subplot(1, 8, 1);
% plot(timeSeries, positionXOf1);
% title('PosX of 1');
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
function a = generateNumericalSolution(starsInSystem, timeSeries, stepSize)
    a = 1;
    % Generate vectors
    positionXOf1 = nan(size(timeSeries));
    positionXOf2 = nan(size(timeSeries));
    positionYOf1 = nan(size(timeSeries));
    positionYOf2 = nan(size(timeSeries));
    
    velocityXOf1 = nan(size(timeSeries));
    velocityXOf2 = nan(size(timeSeries));
    velocityYOf1 = nan(size(timeSeries));
    velocityYOf2 = nan(size(timeSeries));

    % Define initial values
    positionXOf1(1) = 10;
    positionXOf2(1) = 2;
    positionYOf1(1) = 10;
    positionYOf2(1) = 2;
    
    velocityXOf1(1) = 0;
    velocityXOf2(1) = 0;
    velocityYOf1(1) = 0;
    velocityYOf2(1) = 0;
    
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

        forceOn1 = (g_factor ./ (vector_r_length .^ 2)) .* unit_vector_r;
        forceXOn1 = forceOn1(1);
        forceYOn1 = forceOn1(2);

        forceXOn2 = - forceXOn1;
        forceYOn2 = - forceYOn1;

        accelerationXOf1 = forceXOn1 ./ star1.Mass;
        accelerationYOf1 = forceYOn1 ./ star1.Mass;

        accelerationXOf2 = forceXOn2 ./ star2.Mass;
        accelerationYOf2 = forceYOn2 ./ star2.Mass;
        
        velocityXOf1(i + 1) = velocityXOf1(i) + (stepSize .* accelerationXOf1);
        velocityYOf1(i + 1) = velocityYOf1(i) + (stepSize .* accelerationYOf1);

        velocityXOf1(i + 1) = velocityXOf1(i) + (stepSize .* accelerationXOf2);
        velocityYOf2(i + 1) = velocityYOf2(i) + (stepSize .* accelerationYOf2);

        positionXOf1(i + 1) = positionXOf1(i) + (stepSize .* velocityXOf1(i));
        positionYOf1(i + 1) = positionYOf1(i) + (stepSize .* velocityYOf1(i));
        
        positionXOf1(i + 1) = positionXOf2(i) + (stepSize .* velocityXOf2(i));
        positionYOf2(i + 1) = positionYOf2(i) + (stepSize .* velocityYOf2(i));
%       velocity(i + 1) = velocity(i) + (stepSize .* acceleration);
%       position(i + 1) = position(i) + (stepSize * velocity(i));
    end
end

% Returns 1x2 array of distances in (x, y) respectively
function r = calculateDistanceBetweenStars(starA, starB)
    % TODO just reference these directly instead of assigned to temporary
    % variable
    starA_x = starA.Position(1);
    starA_y = starA.Position(2);
    
    starB_x = starB.Position(1);
    starB_y = starB.Position(2);

    r = [starB_x - starA_x, starB_y - starA_y];
end

