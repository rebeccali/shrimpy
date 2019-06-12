function systemStates = myfun(t,zeta)


%air density kg/m^3
rho = 1.225;
%cross sec area m^2
A = 16e-6;
%drag coeff
Cd = 110;
%propeller min radius (m)
Rmin = 0.0035;
%prop max radius
Rmax = 0.0195;
%torque(N*m)
motorTorque = 0.0005;

%Moment of inertia of propeller (Kg*m^2)
Izz = 2*((2.01e-9)+(0.0001*3.5^2)*10^-6);
%Inertial tensor of entire system about COM
inertialT = [ Izz*2000 0 0; 0 Izz*2000 0; 0 0 Izz*1000];
    CdZ = 450;
    mass = 0.004;
    %coefficient of thrust
    kT = 0.114;
    %thrust to COM moment arm distance
    %prop diameter
    propDiam = 0.035;
    %prop Area
    propArea = pi*(propDiam/2)^2;
    %thrust vector needs to update with orientation of the model
    %the thrust vector is always along the local Z axis
    %generate a composite rotation matrix
    %multiply direction vector [0 0 1]' by rotation matrix 
    
    RotY = [cos(zeta(5)), 0, sin(zeta(5));0, 1, 0; -sin(zeta(5)), 0, cos(zeta(5))];
    %Rotation in the Z axis to model effects of precession 
    RotZ = [cos(zeta(5)), -sin(zeta(5)), 0; sin(zeta(5)), cos(zeta(5)),0; 0, 0, 1];
    
    Rcomp= RotY*RotZ; 
    arm = [0.001 0 0];
    arm = Rcomp*arm';
    thrustDirec = Rcomp*[0 0 1]'; 
    thrustDirec = thrustDirec./sqrt(thrustDirec(1)^2+thrustDirec(2)^2+thrustDirec(3)^2);
    %thrustDirec = [0 0 1];
    thrust = (thrustDirec)*kT*rho*(zeta(2).^2)*propDiam^4;
    
    
    %moment about COM due to eccentric thrust
    tMoment = cross([0, 0, norm(thrust)],arm)';  
    
    %calculate angular acceleration of drone due to thrust induced torque
    alpha = inertialT\tMoment; %alpha is a 3x1 vector of angular accelerations about x y z axis
    
    %only consider angular momentum about Z axis
    %propAngMomentum = Izz*zeta(2);
    %propPrecessFreq = tMoment(2)/propAngMomentum;
    %propPrecessFreq = 1;
    %calculate Z axis air drag  
    thrustdragZ = -CdZ*rho*(zeta(4)^2)*propArea/2;
    %thrustdragZ = -CdZ*(zeta(4)^2);
    %calculate drag torque on prop
    dragTq = ((Cd*rho*A/2)*(zeta(2)^2)*(Rmax^3-Rmin^3))/3;
    %dragTq = Cd*(zeta(2)^2);
    %calculate air drag on body
    
    systemStates = [zeta(2); (motorTorque-dragTq)/Izz; zeta(4); (thrust(3) + thrustdragZ - mass*9.81)/mass;zeta(6);alpha(2)];

    %systemStates = [zeta(2); (motorTorque-dragTq)/Izz; zeta(4); (thrust(3) + thrustdragZ - mass*9.81)/mass; 0;0];
    %systemStates = [zeta(2); (motorTorque-dragTq)/Izz; 0;0; 0;0];


end