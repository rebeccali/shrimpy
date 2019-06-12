classdef Propeller < Effector
    %The propeller class remembers everything about a propeller disc
    
    properties
        verbose
        blades
        location
        radius
        beta
        chord
        Cl
        Cd
        axis
    end
    
    methods
        function prop = Propeller(blades, location, radius, beta, chord, axis, varargin)
            try
                prop.verbose = getArgVal('verbose', varargin, false); % Boolean determines if the user receives updates via text on the console
            catch err
                if strcmp(err.identifier,'MATLAB:UndefinedFunction')
                    error('Cannot find getArgVal in the path.  It might be at http://svn.modlabupenn.org/libraries/Matlab_Utilities (rev 67). Find the directory it is in and use addpath() to make it visible.')
                end
            end
            Cl_string = getArgVal('Cl', varargin, 'Cl_segmented'); % The data for the coefficient of lift
            Cd_string = getArgVal('Cd', varargin, 'Cd_segmented'); % The data for the coefficient of drag
            prop.blades = blades;
            prop.location = location;
            prop.radius = radius;
            prop.beta = beta;
            prop.chord = chord;
            prop.axis = axis;
            Cl_in = load(Cl_string);
            prop.Cl = Cl_in.Cl_segmented;
            Cd_in = load(Cd_string);
            prop.Cd = abs(Cd_in.Cd_segmented);
            
        end
        
        function [F, M] = FM(prop, Vcg, Rv_b, omg_cg, omg_bod)
            % Find velocity of center of propeller
%             omg_cg = prop.body.vehicle.state(10:12);
%             omg_bod = prop.body.state(4:6);
%             V = prop.body.vehicle(7:9) + cross3(omg_cg,prop.location);
            V = Vcg + cross3(omg_cg,prop.location);
            % Rotation between vehicle to body
%             Rv_b = RotationMatrix(prop.body.state(1:3));
            F = [0;0;0];
            M = [0;0;0];
            % For each blade
            for blade_num = 1:prop.blades% for each blade
                %Find rotation
                blade_angle = 2*pi*blade_num/prop.blades;
                % Calc rotation matrix between vehicle to propeller
                Rb_p = RotationMatrix(prop.axis*blade_angle);
                Rp_b = Rb_p';
                Rv_p = Rv_b*Rb_p;
                Rp_v = Rv_p';
                % Convert forces to vehicle frame
                [F_temp, M_temp] = prop.ComputeBladeForces(Rp_v*V,Rp_v*omg_cg+Rp_b*omg_bod,0,0,1.225);
                F = F + Rv_p*F_temp;
                M = M + Rv_p*M_temp;
            end
            % Add F*d moments
            M = M + cross3(prop.location,F);
        end
        
        %compute propeller's V and omg from Vcg and flyer omg
        function [F, M] = ComputePropellerForces(prop, angle, omega, V_cg, omg_cg, nu, pitches, rho)
            %Find velocity of the center of the propeller
            V = V_cg + cross3(omg_cg,prop.location);
            F = [0,0,0];
            M = [0,0,0];
            for blade_num = 1:prop.blades% for each blade
                %Find rotation
                blade_angle = 2*pi*blade_num/prop.blades;
                Rr_f = [cos(angle + blade_angle), -sin(angle + blade_angle), 0; sin(angle + blade_angle), cos(angle + blade_angle), 0; 0, 0, 1];
                [Ftemp, Mtemp] = prop.ComputeBladeForces(Rr_f*V,Rr_f*(omg_cg+[0; 0; omega]), nu, pitches(blade_num), rho);
                F = F+Ftemp'*Rr_f;
                M = M+Mtemp'*Rr_f;
            end
            M = M + cross3(prop.location,F);
        end
        
        %Compute the forces and moments of a single blade
        function [F, M] = ComputeBladeForces(prop, V, omg, nu, pitch, rho)
            %ComputeBladeForces(V, omg, nu, pitch, rho) calculates the forces and moments that a
            %blades exerts on the rotor hub.  V is te velocity of the rotor
            %hub in blade coordinates, omg is the angular rate in blade
            %coordinates, nu is the inflow (assumed z), pitch is the blade
            %pitch, and rho is the density
            blade_steps = length(prop.beta);
            F_i = zeros(3,blade_steps);
            M_i = zeros(3,blade_steps);
            step_size = (prop.radius/blade_steps);
            for blade_step = blade_steps:-1:1 %Start from blade tip
                r=blade_step*step_size;
                
                %calculate relative wind
                %TODO:: V is velocity at center of propeller, not this point along blade
%                 Ut1 = V(1) + omg(2)*prop.location(3) - r*omg(3); %correct
                Ut1 = V(1) - r*omg(3);
                Up1 = V(3) + r*omg(1) - nu;
                
                %lift and drag
                phi = atan2(Up1,Ut1);
                aoa1=mod(pi+prop.beta(blade_step)+pitch+phi,2*pi)-pi;
                l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*prop.chord(blade_step)*prop.Cl(int32((aoa1+pi)*100+1));
                d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*prop.chord(blade_step)*prop.Cd(int32((aoa1+pi)*100+1)); % Abs is done at load
                
                %p and t components of force (n is up, t is back)
                n1 = l1*cos(phi)+d1*sin(phi);
                t1 = -l1*sin(phi)+d1*cos(phi);
                
                %force vector
                F_i(1,blade_step) = -t1*step_size;
                F_i(3,blade_step) = -n1*step_size;

                %moment vector
                %Moments from being offset in xyz done in ComputePropellerForces
                M_i(1,blade_step) = - r*n1*step_size;
%                 M_i(2,blade_step) = prop.location(3)*t1*step_size;
                M_i(3,blade_step) = r*t1*step_size;
            end
            F = trapz(F_i,2);
            M = trapz(M_i,2);
%             F = sum(F_i,2);
%             M = sum(M_i,2);
        end
    end
    
end

