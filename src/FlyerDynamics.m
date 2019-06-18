function systemStates = FlyerDynamics(t,zeta)


%air density kg/m^3
rho = 1.225;
%cross sec area m^2
A = 16e-6;
%drag coeff
Cd = 210;
%propeller min radius (m)
Rmin = 0.0035;
%prop max radius
Rmax = 0.0195;
%torque(N*m)
motorTorque = 0.0005;

%Moment of inertia of propeller (Kg*m^2)
Izz = 2*((2.01e-9)+(0.0001*3.5^2)*10^-6);
%Inertial tensor of entire system about COM
inertialT = [ (Izz*70) 0 0; 0 (Izz*70) 0; 0 0 (Izz*10)];
    CdZ = 450;
    %lateral movement drag
    CdL = 2;
    mass = 0.004;
    %coefficient of thrust
    kT = 0.2168;
    %thrust to COM moment arm distance
    %prop diameter
    propDiam = 0.035;
    %prop Area
    propArea = pi*(propDiam/2)^2;
    %thrust vector needs to update with orientation of the model
    %the thrust vector is always along the local Z axis
    %generate a composite rotation matrix
    %multiply direction vector [0 0 1]' by rotation matrix 
    %RotX = [1, 0, 0; 0, cos(zeta(7)), -sin(zeta(7)); 0, sin(zeta(7)),cos(zeta(7))];
    RotY = [cos(zeta(7)), 0, sin(zeta(7));0, 1, 0; -sin(zeta(7)), 0, cos(zeta(7))];
    %Rotation in the Z axis to model effects of precession 
    RotZ = [cos(zeta(5)), -sin(zeta(5)), 0; sin(zeta(5)), cos(zeta(5)),0; 0, 0, 1];
    
    thrustArm = [0.003 0 0]'; % eccentric thrust moment arm
    diffLiftArm = [0.005 0 0]'; % differential lift moment arm
    thrustDirec = (RotZ*RotY*[0 0 1]');
    thrust = (thrustDirec)*(kT*rho*(zeta(2)^2)*propDiam^4);

    diffLift = (2*(zeta(7)^2)*diffLiftArm(1));

    %moment about COM due to eccentric thrust
    tMoment = RotZ*((cross(thrustArm, thrust))-diffLift);
    %tMoment = [tMoment(1),tMoment(2)+diffLift,tMoment(3)]';
    %calculate angular acceleration of drone due to thrust induced torque
 
    %only consider angular momentum about Z axis
    
    
    %calculate Z axis air drag  
    thrustdragZ = -CdZ*rho*(zeta(4)^2)*propArea/2;
    %thrustdragZ = -CdZ*(zeta(4)^2);
    %calculate drag torque on prop
    dragTq = ((Cd*rho*A/2)*(zeta(2)^2)*(Rmax^3-Rmin^3))/3;
    alphaBody = inertialT\[0 0 (-motorTorque+dragTq)]';
    alphaPitch = (inertialT\(tMoment));
    
    %lateral thrusts in the X and Y are not equal nor does it cancel out
    %through rotation. Not sure why. :(((( Has to do with effects of
    %differential lift 
    latThrustX = (thrust(1) -(CdL*zeta(10)))/mass;
    latThrustY = (thrust(2) - (CdL*zeta(12)))/mass;

    
    %alphaBody = (-motorTorque+dragTq)/(Izz*10);
    
    %dragTq = Cd*(zeta(2)^2);
    %calculate air drag on body
    %The states are organized as: 
    %[propeller Angular position; propeller angular velocity; Z axis position; Z axis velocity;...
    % body angular position; body angular velocity; pitch angular position, pitch angular velocity;...
    % X positon, X velocity, Y position, Y velocity]
    systemStates = [zeta(2); (motorTorque-dragTq)/Izz; zeta(4); (thrust(3) + thrustdragZ - mass*9.81)/mass;...
        zeta(6); alphaBody(3); zeta(8); alphaPitch(2); zeta(10); latThrustX; zeta(12); latThrustY];

    %systemStates = [zeta(2); (motorTorque-dragTq)/Izz; zeta(4); (thrust(3) + thrustdragZ - mass*9.81)/mass; 0;0];
    %systemStates = [zeta(2); (motorTorque-dragTq)/Izz; 0;0; 0;0];


end