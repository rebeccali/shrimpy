classdef VehicleElement < handle
    %A VehicleElement is an object that the vehicle class deals with
    %directly.
    
    properties
        state
    end
    
    methods (Abstract)
        [dX, F, M, I] = ODEFMI(obj,V_cg,omg_cg); % ODE with net (external) forces and moments
    end
    
end

