classdef Vehicle < handle
    %The Vehicle class describes a specific physical flying vehicle
    
    properties
        verbose
        elements
        m
        I_static
        state % [x, y, z, phi, theta, psi, x_dot, y_dot, z_dot, p, q, r] ... 
                % x y z phi theta psi are in world coordinates
                % x_dot y_dot z_dot p q r are in vehicle coordinates
    end
    
    properties (Access = private)
        state_elements
    end
    
    methods
        function veh = Vehicle(varargin)
            %Vehicle(bodies, varargin) is a constructor for a vehicle.  
            %bodies is an array of bodies that comprise the vehicle
            try
                veh.verbose = getArgVal('verbose', varargin, false); % Boolean determines if the user receives updates via text on the console
            catch err
                if strcmp(err.identifier,'MATLAB:UndefinedFunction')
                    error('Cannot find getArgVal in the path.  It might be at http://svn.modlabupenn.org/libraries/Matlab_Utilities (rev 67). Find the directory it is in and use addpath() to make it visible.')
                end
            end
            veh.state = getArgVal('state', varargin, [0,0,0,0,0,0,0,0,0,0,0,0]); % The initial state of the vehicle
            veh.state_elements = 12;
            veh.elements = [];
            veh.m = 0;
            veh.I_static = zeros(3);
            
            elements = getArgVal('elements', varargin, []); % In case you want to pass in the elements
            if(~isempty(elements))
                RegisterElements(veh,elements);
            end
        end
        
        function RegisterElements(veh, elements)
            veh.elements = [veh.elements elements];
            veh.state_elements = 12; % This vehicle is also sort of a vehicle elements in that it has 12 states
            for i = 1:length(veh.elements)
                veh.state_elements(i+1) = length(veh.elements(i).state);
            end
            for i = isa(veh.elements, 'Body')
                veh.m = veh.m + veh.elements(i).m;
                veh.I_static = veh.I_static + veh.elements(i).m*((veh.elements(i).location')*veh.elements(i).location*eye(3)-veh.elements(i).location*(veh.elements(i).location'));
            end
        end
        
        function dX = ODE(veh,t,X)
            %% Update motor and body states
            veh.state = X(1:12);
%             pos = X(1:3);
            ang = X(4:6);
            V_cg = X(7:9);
            omg_cg = X(10:12);
            dX = zeros(sum(veh.state_elements),1);
            element_count = 12;
            for i = 1:length(veh.elements)
                veh.elements(i).state = X(element_count+1:element_count+veh.state_elements(i+1));
                element_count = element_count + veh.state_elements(i);
            end
            %% Cycle through elements (motors first by way elements packed into veh.elements
            F = [0;0;0];
            M = [0;0;0];
            I = veh.I_static;
            element_count = 12;
            for i = 1:length(veh.elements)
                % Compute ODE forces and moments
                [dX_temp, F_temp, M_temp, I_temp] = veh.elements(i).ODEFMI(V_cg,omg_cg);
                % Add the resulting diff to dX
                dX(element_count+1:element_count+veh.state_elements(i+1)) = dX_temp;
                % Add the resulting F and M to the net F and M
                F = F + F_temp;
                M = M + M_temp;
                I = I + I_temp;
            end
            %% Run vehicle ODE
            % Calculate rotation matrix
            Rv_i = RotationMatrix(ang)'; % ang is world to vehicle, so transpose returned Ri_v to get Rv_i.
            % Calculate CG acceleration
            dX(7:9) = F/veh.m-cross3(omg_cg,V_cg) + Rv_i*[0;0;9.81];
            % Calculate CG angular acceleration
            dX(10:12) = I\M;
            % Calculate World angular velocity (this is just transformed vehicle rates
            dX(4:6) = BodyToWorldRates(omg_cg,ang);
            % Calculate World Velocity (this is just transformed body velocity
            dX(1:3) = Rv_i'*V_cg; % Vcg is in vehicle frame
            % Assemble state_dot
        end
        
        function stop = Control(veh, t, X, flag) 
            if strcmp(flag,'init')
                
            elseif strcmp(flag,'done')

            else
                stop = false;
                % Compute motor voltages
            end
        end
    end
    
end

