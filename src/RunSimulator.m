%% Runs a simulation of the vehicle and body dynamics 

%% Clean Space
close all;
clear all;

%% Define Simulation Parameters
dt = 0.01;
startTime = 0;
endTime = 20;
x0 = [0;0;0;0.5;0;200;0;0;0;0;0;0];

%% Do Simulation

[t,states] = ode45(@FlyerDynamics, [startTime, endTime], x0);

%% Plot states
% plotFlyer(t, states);

%% Do Animation
animateFlyer(t, states, dt, startTime, endTime)











