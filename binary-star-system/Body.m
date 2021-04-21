classdef Body
    %   Body
    %   This class represents a cellestial body. A body can be a star, a
    %   small moon, or a grain of sand.
    
    % Properties that can only be modified by the class constructor
    properties (SetAccess = immutable)
        Weight double {mustBeReal, mustBeFinite, mustBeNonnegative} % in kg
    end

    properties
        % (1, 1) is x coordinate
        % (1, 2) is y coordinate
        Position(1,2) double {mustBeReal, mustBeFinite} % in metres?
        Velocity(1,2) double {mustBeReal, mustBeFinite} % in metres?
    end

    methods
        % Constructor of Body class
        function obj = Body(weight,position, velocity)
            obj.Weight = weight;
            obj.Position = position;
            obj.Velocity = velocity;
        end
    end
end

