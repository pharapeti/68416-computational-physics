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
        Bodies
    end
    
    methods
        function obj = Simulation(time_step,total_time, G)
            obj.TimeStep = time_step;
            obj.TotalTime = total_time;
            obj.TimeSeries = 0:time_step:total_time;
            obj.G = G;
        end
       
        function obj = addBody(obj, body)
            assert(isa(body,'Body'),...
                'Individual Constructor Error:  body is of class %s, not a Body object.', class(body));

            obj.Bodies = [obj.Bodies, body];
        end
    end
end

