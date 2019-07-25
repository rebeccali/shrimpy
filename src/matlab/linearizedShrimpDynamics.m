function dxdt = linearizedShrimpDynamics(x, t)
%% Linearized Shrimp lateral dynamics based off of the model in Piccoli's Thesis
%  Reference page 19, equation 2.11. 
%  Inputs:
%  x : state vector [u, v, p q, phi, theta]
%  t : time (unused)
%  Outputs:
%  dxdt : derivative of state for ODE integration 
%
%  State definitions:
%  Flyer frame is the frame where the yaw about the axis of rotation is 
%  fixed.
%  u = x-velocity of the flyer in the inertial frame (p.15)
%  v = y-velocity of the flyer in the inertial frame (p.15)
%  p = x-axis angular velocity of flyer in the inertial frame (p.15)
%  q = y-axis angular velocity of flyer in the inertial frame (p.15)
%  phi = x-axis angle world frame to flyer frame (pitch)
%  theta = y-axis angle world frame to flyer frame (roll)
%  Forces and Moments:
%  f_a = [X;Y;Z]; which is forces in flyer frame (p.15)
%  tau_a = [L M N]'; which is moment in flyer frame (P.15)

u = x(1);
v = x(2);
p = x(3);
q = x(4);
phi = x(5);
theta = x(6);

v_flyer_w = [u v 0]'; % xyz velocity of flyer in inertial frame 
r_differentialLift_f = ???; % Differential lift in the flyer frame 

% Definte parts of the system and environment
Ixx = 1;
Iyy = 1;
Izz = 2;
mass = 1;
g = 9.81; % Gravitational Constant [m/s^2]

% Define our aerodynamics
copComMomentPerpToV = 0; % COP = COM
gyroPrecessionMomentPerpToOmega = = -Izz*r/Ixx;



% Definte all forces and moments.
% V is the lateral velocity vector
% Omega is the angular velocity of the shaft rotation vector
a = dragForceParallelToV;
b = diffLiftMomentParallelToV;
c = copComMomentPerpToV;
d = dragMomentParallelToOmega;
e = gyroPrecessionMomentPerpToOmega;

% Assertions from p.19
assert(a<0.0, "a should always be negative");
assert(d<0.0, "d should always be negative");



% Jacobian/State transition matrix 
jacobian = [ a  0  0  0  0 -g;
             0  a  0  0  g  0;
             b  c  d  e  0  0;
             -c b -e  d  0  0;
             0  0  1  0  0  0;
             0  0  0  1  0  0;]

end 