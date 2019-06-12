classdef Body < VehicleElement
    %The Body class describes a rigid body.  A body has a mass, inertia, a
    %location (relative to the vehicle CG), and optionally a propeller
    
    properties
        verbose
        m
        I
        location
        effectors
        link % each body constrained to only be connected to one link
        % state % [phi, theta, psi, p, q, r] % phi, theta, psi in Vehicle frame, p, q, r in body frame
    end
    
    methods
        function bod = Body(mass, inertia, location, effectors, link, varargin)
            try
                bod.verbose = getArgVal('verbose', varargin, false); % Boolean determines if the user receives updates via text on the console
            catch err
                if strcmp(err.identifier,'MATLAB:UndefinedFunction')
                    error('Cannot find getArgVal in the path.  It might be at http://svn.modlabupenn.org/libraries/Matlab_Utilities (rev 67). Find the directory it is in and use addpath() to make it visible.')
                end
            end
            bod.state = getArgVal('state', varargin, [0,0,0,0,0,0]); % The initial state of the body
            bod.m = mass;
            bod.I = inertia;     
            bod.location = location;
            bod.effectors = effectors;
            bod.link = link;
        end
        
        function [F, M, I] = FMI(bod,V_cg,omg_cg,ang_bod,omg_bod)
            F = [0;0;0];
            M = [0;0;0];
            % Calculate rotation matrix between body and vehicle
            Rv_b = RotationMatrix(ang_bod); % phi, theta, psi from vehicle frame to body frame, gives Rv_b.
            % Compute forces from effectors
            for i = 1:length(bod.effectors)
                [F_temp, M_temp] = bod.effectors(i).FM(V_cg, Rv_b, omg_cg, omg_bod);
                F = F+F_temp;
                M = M+M_temp;
            end
            % Calculate gyroscopic forces
            M = M - cross3(omg_cg,bod.I*(omg_cg+omg_bod)); % Not missing I*omg_dot, that is for the "untransferred moments"
            % Calculate instantaneous vehicle inertia
            I = Rv_b*bod.I*Rv_b'; % TODO:: CHECK MY MATH!!!!!!!!!!
        end
%         
%         function [dX, F, M, I] = ODEFMI(bod,V_cg,omg_cg)
%             F_gen = [0;0;0];
%             M_gen = [0;0;0];
%             % Calculate rotation matrix between body and vehicle
%             Rv_b = RotationMatrix(bod.state(1:3)); % phi, theta, psi from vehicle frame to body frame, gives Rv_b.
%             % Compute forces from effectors
%             for i = 1:length(bod.effectors)
%                 [F_temp, M_temp] = bod.effectors(i).FM(V_cg, Rv_b, omg_cg, bod.state(4:6));
%                 F_gen = F_gen+F_temp;
%                 M_gen = M_gen+M_temp;
%             end
%             % Calculate gyroscopic forces
%             M_gen = M_gen - cross3(omg_cg,bod.I*(omg_cg+bod.state(4:6))); % are these omgs in the same frame??? % Not missing I*omg_dot, that is for the "untransferred moments"
%             % Transfer forces through links
%             [F, M] = bod.link.TransferredFM(F_gen,M_gen);
%             % Subtract off transferred M to find net moments for this body
%             M_body = M_gen - M;
%             % Calculate dX for p, q, r
%             % Calculate phi_dot, theta_dot, psi_dot, which are rotated p, q, r
%             dX = [BodyToWorldRates(bod.state(4:6),bod.state(1:3));bod.I\Rv_b'*M_body];
%             % Calculate instantaneous vehicle inertia
%             I = Rv_b*bod.I*Rv_b'; % TODO:: CHECK MY MATH!!!!!!!!!!
%         end
    end
    
end

