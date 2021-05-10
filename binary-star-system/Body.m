classdef Body
    %   Body
    %   This class represents a cellestial body. A body can be a star, a
    %   small moon, or a grain of sand.
    
    % Properties that can only be modified by the class constructor
    properties (SetAccess = immutable)
        Mass double {mustBeReal, mustBeFinite, mustBeNonnegative} % in kg
    end

    properties
        Position % m
        Velocity % m/s
        Acceleration % m/s/s
        KE
        PE
    end

    methods
        % Constructor of Body class
        function obj = Body(mass, initial_position, initial_velocity, initial_acceleration, time_series)
            obj.Mass = mass;

            % Use time series to generate position and velocity vectors for
            % each time iteration
            obj.Position.X = nan(1, length(time_series));
            obj.Position.Y = nan(1, length(time_series));

            obj.Velocity.X = nan(1, length(time_series));
            obj.Velocity.Y = nan(1, length(time_series));
            
            obj.Acceleration.X = nan(1, length(time_series));
            obj.Acceleration.Y = nan(1, length(time_series));

            obj.KE = nan(1, length(time_series));
            obj.PE = nan(1, length(time_series));

            % Set initial position and velocity values
            obj.Position.X(1) = initial_position(1);
            obj.Position.Y(1) = initial_position(2);

            obj.Velocity.X(1) = initial_velocity(1);
            obj.Velocity.Y(1) = initial_velocity(2);
            
            obj.Acceleration.X(1) = initial_acceleration(1);
            obj.Acceleration.Y(1) = initial_acceleration(2);

            % TODO: Calculate this based on initial conditions
            obj.KE(1) = 0;
            obj.PE(1) = 0;
        end
    end
end

