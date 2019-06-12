function [tout, Xout] = LinearizedSimulator( X, time, A )
%LinearizedSimulator Simulates the linearized version of yimFlyerLite
%Syntax: [tout, Xout] = LinearizedSimulator( X, time, A )
% Where dot(X) = AX

[tout, Xout] = ode113(@LinearizedSimulatorODE,linspace(0,time,time/.001),X,odeset('MaxStep',.001,'OutputFcn',@controlLoop));

function stop = controlLoop(t, X, flag) 
        if strcmp(flag,'init')

        elseif strcmp(flag,'done')

        else
            
        end
end

function dX = LinearizedSimulatorODE(t,X)
    dX = A*X;
end
end

