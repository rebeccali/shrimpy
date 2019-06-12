classdef RigidLink < Link
    %The Motor class describes and simulates a motor
    
    properties
        verbose
    end
    
    methods
        function lnk = RigidLink(varargin)
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
            lnk.state = [];
        end
        
        function [dX, F, M, I] = ODEFMI(lnk)
            dX = [];
            F = [0, 0, 0];
            M = [0, 0, 0];
            I = zeros(3);
        end
        
        function [F_out, M_out] = TransferredFM(lnk, F_in, M_in)
            F_out = F_in;
            M_out = M_in;
        end
    end
    
end

