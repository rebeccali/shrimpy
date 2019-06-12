% Compares nonlinear sim to linear sim.  Non-linear sim must be included in
% path.  Compares pulsing.
clear all; close all; clc;
sim_time = 10;
Piccolissimo_V11;
disturbances = [nan; .5; .5; nan; .05; .05; nan; nan; nan; .1; .1; nan; nan; nan; nan; nan; nan; nan];
[state_transition, inds_performed] = GenerateStateTransition(disturbances);
motor_offset = .00807; args = {'throttle', 'cyclic', 'plot'};
[tout, Xout] = FlyerSimulation(Xbase,sim_time,args);
[t,X] = simulatePiccolissimoV1(state_transition,'plot_hold',true,'t_max',sim_time);
beep;