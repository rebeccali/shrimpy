%% Run a simple simulation for Piccolissimo that works.
%  This script is used fo regression testing, make a separate
%  script for experimentation.
%  Copyright Modlab 2019
%  Author: Rebecca Li

close all;
clear all;
% Add all folders
addpath(genpath('./'));

% Initialize global parameters Piccolissimo
Piccolissimo_V11;

% Run a test simulation
global Xbase
X = Xbase;
startTime = 0;
finalTime = 0.5;
time_in = [startTime, finalTime];
argsIn = {'throttle', 'circle','plot'};
[tout, Xout] = FlyerSimulation(X,time_in,argsIn);

% Write simulation to test files. If they have changed, we 
% have broken the simulation somehow.
csvwrite('test/testXout.csv', Xout);
csvwrite('test/testTout.csv', tout);