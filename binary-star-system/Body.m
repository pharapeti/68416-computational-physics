classdef Body
    %   Body
    %   This class represents a cellestial body. A body can be a star, a
    %   small moon, or a grain of sand.
    
    % Properties that can only be modified by the class constructor
    properties (SetAccess = immutable)
        Mass double {mustBeReal, mustBeFinite, mustBeNonnegative} % in kg
    end

    properties
        Position double % in metres?
        Velocity double % in metres?
    end

    methods
        % Constructor of Body class
        function obj = Body(mass, initial_position, initial_velocity, time_series)
            obj.Mass = mass;

            % Use time series to generate position and velocity vectors for
            % each time iteration
            obj.Position = nan(length(time_series), 2);
            obj.Velocity = nan(length(time_series), 2);

            % Set initial position and velocity values
            obj.Position(1, :) = initial_position;
            obj.Velocity(1, :) = initial_velocity;
        end
    end
end

