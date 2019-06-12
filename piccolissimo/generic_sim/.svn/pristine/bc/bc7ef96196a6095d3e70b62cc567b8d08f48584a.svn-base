%Test vehicle setup
% Don't trust this script to be anything, it changes
% clear all;
%% Make the vehicle
veh = Vehicle();

%% Make the propeller
prop = Propeller(2, [0;0;0], .140, atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]), [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049], [0,0,1]);

%% Make the links
lnk = RigidLink(veh);

%% Make the bodies
bod = Body(.039, [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9, [0;0;0], prop, lnk);

%% Register vehicle elements
veh.RegisterElements(bod);

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
