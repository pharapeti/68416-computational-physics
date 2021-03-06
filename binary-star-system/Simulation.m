classdef Simulation < handle
    %   Simulation
    %   The Simulation class defines the system parameters relating the
    %   physical aspects of the system, including the simulation
    %   characteristics and values such as the Gravitational Constant.

    % The following properties can only be set by the constuctor
    properties (SetAccess = immutable)
        TimeStep double {mustBeReal, mustBeFinite, mustBeNonnegative} % seconds
        TotalTime int32 {mustBeReal, mustBeFinite, mustBeNonnegative} % seconds
        TimeSeries(1,:) {mustBeReal, mustBeFinite, mustBeNonnegative} % seconds
        G double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    properties
        Bodies        % array of Bodies in the simulation
        Barycenter    % (x, y) for in iteration of the simulation
    end

    methods
        function obj = Simulation(time_step,total_time, G)
            obj.TimeStep = time_step;
            obj.TotalTime = total_time;
            obj.TimeSeries = 0:time_step:total_time;
            obj.G = G;

            obj.Bodies = [];
            obj.Barycenter = Barycenter();
            obj.Barycenter.X = nan(1, length(obj.TimeSeries));
            obj.Barycenter.Y = nan(1, length(obj.TimeSeries));
        end

        function obj = createBody(obj, mass, init_position, init_velocity, init_acceleration)
            newBody = Body(mass, init_position, init_velocity, init_acceleration, obj.TimeSeries);
            obj.Bodies = [obj.Bodies, newBody];
        end
    end
end

