clc; clear; close all;

initial_height = 1000; % metres
gravity = -9.8; % metres/second/second
dt = 0.01; % seconds
timeSeries = 0:dt:20;

height = nan(size(timeSeries));
height(1) = initial_height;

velocity = nan(size(timeSeries));
velocity(1) = 0;

% velocity_new = velocity_old + (acceleration * time);
% height_new = height_old + (dt * velocity);

time = nan; % when ball hits ground

for n = 1:length(height)-1    
    velocity(n+1) = velocity(n) + (gravity .* dt);
    height(n+1) = height(n) + (dt * velocity(n));
    
    if height(n+1) <= 0
        time = timeSeries(n);
        break;
    end
    
end

disp(time);
plot(timeSeries, height);