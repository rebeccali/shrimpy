%% Run a simple simulation for Piccolissimo that works.
%  Copyright Modlab 2019
%  Author: Rebecca Li

close all;
clear all;
% Add all folders
addpath(genpath('./'));

% Initialize global parameters Piccolissimo
Piccolissimo_V11;

% Run a simulation
global Xbase
X = Xbase;
startTime = 0;
finalTime = 10;
time_in = [startTime, finalTime];
argsIn = {'throttle', 'cyclic','plot'};
[tout, Xout] = FlyerSimulation(X,time_in,argsIn)