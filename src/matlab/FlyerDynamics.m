function dzetadt = AllenFlyerDynamics(t, zeta)


%air density kg/m^3
rho = 1.225;
%cross sec area m^2
surfAreaProp = 16e-6;
%drag coeff
coeffDragProp = 210; % TODO replace
%propeller min radius (m)
radiusPropMin = 0.0035;
%prop max radius
radiusPropMax = 0.0195;
%torque(N*m)
motorTorque = 0.0005;

%Moment of inertia of propeller (Kg*m^2)
Izz = 2*((2.01e-9)+(0.0001*3.5^2)*10^-6);
%Inertial tensor of entire system about COM
inertialT = [ (Izz*70) 0 0; 0 (Izz*70) 0; 0 0 (Izz*10)]; % ?? entire system or system minus prop?
    coeffDragBodyZ = 450; % 
    %lateral movement drag
    coeffDragBodyXY = 2;
    mass = 0.004;
    %coefficient of thrust
    kT = 0.2168;
    %thrust to COM moment arm distance
    %prop diameter
    propDiam = 0.035;
    %prop Area
    propDiscArea = pi*(propDiam/2)^2;
    %thrust vector needs to update with orientation of the model
    %the thrust vector is always along the local Z axis
    %generate a composite rotation matrix
    
    %multiply direction vector [0 0 1]' by rotation matrix 
    %RotX TODO
    RotY = [cos(zeta(7)), 0, sin(zeta(7));0, 1, 0; -sin(zeta(7)), 0, cos(zeta(7))];
    %Rotation in the Z axis to model effects of precession 
    RotZ = [cos(zeta(5)), -sin(zeta(5)), 0; sin(zeta(5)), cos(zeta(5)),0; 0, 0, 1];
    rotmRotationAxis2World = RotZ*RotY;
    
    thrustArm2Body_Body = [0.003 0 0]'; % eccentric thrust moment arm
    thrustArm2Body_World = TODO;
    diffLiftArm2Body_Body = [0.005 0 0]'; % differential lift moment arm
    thrustDirection_Body = [0 0 1]';
    thrustDirection_World = rotmRotationAxis2World * thrustDirection_Body;
    thrustForceVector_World = thrustDirection_World*(kT*rho*(zeta(2)^2)*propDiam^4); % online 

    diffLiftMoment_World = (2*(zeta(7)^2)*diffLiftArm2Body_Body(1)); % TODO replace with piccoli calc

    %moment about COM due to eccentric thrust
    thrustMoment_World = cross(thrustArm2Body_World, thrustForceVector_World)-diffLiftMoment_World;
    %calculate angular acceleration of drone due to thrust induced torque
 
    
    %calculate Z axis air drag  
    thrustdragZ = -coeffDragBodyZ*rho*(zeta(4)^2)*propDiscArea/2;
    
    %calculate drag torque on prop
    dragTq = ((coeffDragProp*rho*surfAreaProp/2)*(zeta(2)^2)*(radiusPropMax^3-radiusPropMin^3))/3;
    angularAccelBodyTensor = inertialT\[0 0 (-motorTorque+dragTq)]'; % TODO
    angularAccelTensor = (inertialT\(thrustMoment_World)); % TODO
    angularAccelPropZ = (motorTorque-dragTq)/Izz; % TODO?
    
    %lateral thrusts in the X and Y are not equal nor does it cancel out
    %through rotation. Not sure why. :(((( Has to do with effects of
    %differential lift 
    accelX = (thrustForceVector_World(1) -(coeffDragBodyXY*zeta(10)))/mass;
    accelY = (thrustForceVector_World(2) - (coeffDragBodyXY*zeta(12)))/mass;
    accelZ = (thrustForceVector_World(3) + thrustdragZ - mass*9.81)/mass;

    
    %calculate air drag on body
    %The states are organized as: 
    %[propeller Angular position; propeller angular velocity; Z axis position; Z axis velocity;...
    % body angular position; body angular velocity; pitch angular position, pitch angular velocity;...
    % X positon, X velocity, Y position, Y velocity]
    dzetadt = [zeta(2); angularAccelPropZ; zeta(4); accelZ;...
        zeta(6); angularAccelBodyTensor(3); zeta(8); angularAccelTensor(2); zeta(10); accelX; zeta(12); accelY];

end