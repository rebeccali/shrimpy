% Quadrotor test setup
% Don't trust this script to be anything, it changes
% clear all;

%% Make the propeller
prop1 = Propeller(2, [.25;0;0], .140, atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]), [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049], [0,0,1]);
prop2 = Propeller(2, [0;-.25;0], .140, atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]), [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049], [0,0,1]);
prop3 = Propeller(2, [-.25;0;0], .140, atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]), [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049], [0,0,1]);
prop4 = Propeller(2, [0;.25;0], .140, atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]), [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049], [0,0,1]);


%% Make the bodies
bod = Body(
bod1 = Body(.039, [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9, [.25;0;0], prop1);
bod2 = Body(.039, [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9, [0;-.25;0], prop2);
bod3 = Body(.039, [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9, [-.25;0;0], prop3);
bod4 = Body(.039, [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9, [0;.25;0], prop4);

%% Make the links
mot1 = Motor(1/(740*2*pi/60),.26,.02,bod1, [0 0 1]);
mot2 = Motor(1/(740*2*pi/60),.26,.02,bod2, [0 0 1]);
mot3 = Motor(1/(740*2*pi/60),.26,.02,bod3, [0 0 1]);
mot4 = Motor(1/(740*2*pi/60),.26,.02,bod4, [0 0 1]);

%% Make the vehicle
veh = Vehicle();
%% Register vehicle elements
veh.RegisterElements([bod1]);

%% Setup and run ODE
time = 3;
initial_vehicle_state = [ ...
    0; ... %vehicle x
    0; ... %vehicle y
    0; ... %vehicle z
    0; ... %vehicle phi
    0; ... %vehicle theta
    0; ... %vehicle psi
    0; ... %vehicle x_dot
    0; ... %vehicle y_dot
    0; ... %vehicle z_dot
    0; ... %vehicle p
    0; ... %vehicle q
    0; ... %vehicle r
    ];

initial_body_state = [ ...
    0; ... %vehicle phi
    0; ... %vehicle theta
    0; ... %vehicle psi
    0; ... %vehicle p
    0; ... %vehicle q
    0; ... %vehicle r
    ];
initial_state = [initial_vehicle_state; initial_body_state];
% dX = veh.ODE(0,initial_state)
[tout, Xout] = ode113(@veh.ODE,linspace(0,time,time/.001),initial_state,odeset('MaxStep',.001)); %,'OutputFcn',@controlLoop
