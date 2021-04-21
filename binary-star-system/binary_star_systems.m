%% Part Zero : Setup

% Clear existing workspace
clear; clc; close all

% Setup parameters
timestep = 0.01; % timestep (seconds)
totalTime = 10; % total time of simulation (seconds)
timeSeries = 0:timestep:totalTime;

% Define initial conditions of system
xInitial = 2; % initial position (metres)
vInitial = 0; % intial velocity (metres / second)

% Constants
G = 6.67 * 10.^-11;
M1 = 10; % kg
M2 = 20; % kg

star1 = Body(10, [5, 5], [0, 0]);
star2 = Body(20, [0, 0], [0, 0]);

[positionXOf1, positionXOf2, positionYOf1, positionYOf2, velocityXOf1, velocityXOf2, velocityYOf1, velocityYOf2] = generateNumericalSolution(timeSeries, timestep);

subplot(1, 8, 1);
plot(timeSeries, positionXOf1);
title('PosX of 1');

subplot(1, 8, 2);
plot(timeSeries, positionXOf2);
title('PosX of 2');

subplot(1, 8, 3);
plot(timeSeries, positionYOf1);
title('PosY of 1');

subplot(1, 8, 4);
plot(timeSeries, positionYOf2);
title('PosY of 2');

subplot(1, 8, 5);
plot(timeSeries, velocityXOf1);
title('VelX of 1');

subplot(1, 8, 6);
plot(timeSeries, velocityXOf2);
title('VelX of 2');

subplot(1, 8, 7);
plot(timeSeries, velocityYOf1);
title('VelY of 1');

subplot(1, 8, 8);
plot(timeSeries, velocityYOf2);
title('VelY of 2');

% This function computes a numerical solution for a given dataset via the
% Euler's method.
function [positionXOf1, positionXOf2, positionYOf1, positionYOf2, velocityXOf1, velocityXOf2, velocityYOf1, velocityYOf2] = generateNumericalSolution(timeSeries, stepSize)
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
    
    % Gravitational Factor
    G = 6.67408 * 10.^-11;
    M1 = 10; % kg
    M2 = 20; % kg
    g_factor = G * M1 * M2;

    % Inspired by http://www.astro.yale.edu/coppi/astro520/solving_differential_equation.pdf
    % Numerically solve the position and velocity of the model system using
    % function arguments
    for i = 1:length(timeSeries) - 1
        % Calculate distance between star 1 and 2
        vector_r = calculateDistanceBetweenStars(3, 2, 1, 5);
        vector_r_length = norm(vector_r);
        unit_vector_r = vector_r ./ vector_r_length;

        forceOn1 = (g_factor ./ (vector_r_length .^ 2)) .* unit_vector_r;
        forceXOn1 = forceOn1(1);
        forceYOn1 = forceOn1(2);

        forceXOn2 = - forceXOn1;
        forceYOn2 = - forceYOn1;

        accelerationXOf1 = forceXOn1 ./ M1;
        accelerationYOf1 = forceYOn1 ./ M1;

        accelerationXOf2 = forceXOn2 ./ M2;
        accelerationYOf2 = forceYOn2 ./ M2;
        
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

function r = calculateDistanceBetweenStars(xpos, ypos, x2pos, y2pos)
    r = [x2pos - xpos, y2pos - ypos];
end