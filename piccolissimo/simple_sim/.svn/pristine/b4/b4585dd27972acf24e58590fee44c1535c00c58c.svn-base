% Compares nonlinear sim to linear sim.  Non-linear sim must be included in
% path.  Compares response from disturbance
clear all; close all; clc;
velocity_disturbance = .1;
sim_time = 10;
Piccolissimo_V11;
disturbances = [nan; velocity_disturbance; velocity_disturbance; nan; .05; .05; nan; nan; nan; .01; .01; nan; nan; nan; nan; nan; nan; nan];
[state_transition, inds_performed] = GenerateStateTransition(disturbances);
motor_offset = .00807; args = {'throttle', 'plot'};
x0 = Xbase + [0; velocity_disturbance; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
[tout, Xout] = FlyerSimulation(x0,sim_time,args);
[t,X] = simulatePiccolissimoV1(state_transition,'plot_hold',true,'t_max',sim_time,'pulsing',false,'x0',[velocity_disturbance,0,0,0,0,0,0,0,0,0,0,0]);
beep;
% Figures 1, 2, 33, 36