classdef Motor < Link
    %The Motor class describes and simulates a motor
    
    properties
        verbose
        Kt
        R
        L
        V
        stator % Body that the stator is attached to
        rotor % Body that the rotor is attached to
        axis
        %state = i
        M % Internal moment (along axis) that is created by the motor
    end
    
    methods
        function mot = Motor(Kt, R, L, stator, rotor, axis, varargin)
            %Motor(Kt, R, L, stator, rotor, axis, varargin) is the motor constructor.
            % Axis is the SO3 axis of rotation (ex. [0 0 1];)
            % Stator and rotor are body objects.
            try
                veh.verbose = getArgVal('verbose', varargin, false); % Boolean determines if the user receives updates via text on the console
            catch err
                if strcmp(err.identifier,'MATLAB:UndefinedFunction')
                    error('Cannot find getArgVal in the path.  It might be at http://svn.modlabupenn.org/libraries/Matlab_Utilities (rev 67). Find the directory it is in and use addpath() to make it visible.')
                end
            end
            mot.state = getArgVal('state', varargin, 0); % The initial state of the motor
            mot.axis = axis/norm(axis);
            mot.Kt = Kt;
            mot.R = R;
            mot.L = L;
            if isa(stator,'Body') && isa(rotor,'Body')
                mot.stator = stator;
                mot.rotor = rotor;
            else
                error('The stator and/or rotor are not Body objects');
            end
        end
        
        function [dX, F, M, I] = ODEFMI(mot)
            dX = (mot.V - mot.Kt*(mot.rotor.state(12)-mot.stator.state(12))-mot.R*mot.state)/mot.L;
            F = [0, 0, 0]; % Motors make inertnal forces. F and M are external forces.
            M = [0, 0, 0];
            mot.M = mot.Kt*X;
            I = zeros(3);
        end
        
        function [F_out, M_out] = TransferredFM(mot, F_in, M_in)
            F_out = F_in;
            % TODO:: CHECK THIS MATH
            M_out = ([1,1,1] - mot.axis).*M_in + mot.M.*mot.axis;
        end
    end
    
end

