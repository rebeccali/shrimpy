clear all
close all
clc

syms xy z
syms w1 w2 w3;
syms M1 M2 M3;

%% Differential Lift
disp('DIFF LIFT!!!!!!!!!!!!!!!!!!!!!')
%% Step One: Advancing/Retreating effect causes moment
disp('Advancing/Retreating effect causes moment');
I = [xy,0,0;0,xy,0;0,0,z];
w = [0;0;w3];
M = [M1;0;0];
wDot = I\(M-cross(w,I*w))

%% Step Two: Gyro term causes moment
disp('Gyro term causes moment');
w = [w1;0;w3];
wDot = I\(M-cross(w,I*w))

%% Step Three: Motion caused by gyro causes new motion by gyro
disp('Motion caused by gyro causes new motion by gyro');
w = [w1;w2;w3];
wDot = I\(M-cross(w,I*w))

%% Drag and CG misaligned
disp('DRAG!!!!!!!!!!!!!!!!!!!!!!!!!!!')
%% Step One: Drag creates moment
disp('Drag creates moment');
I = [xy,0,0;0,xy,0;0,0,z];
w = [0;0;w3];
M = [0;M2;0];
wDot = I\(M-cross(w,I*w))

%% Step Two: Gyro term causes moment
disp('Gyro term causes moment');
w = [0;w2;w3];
wDot = I\(M-cross(w,I*w))

%% Step Three: Motion caused by gyro causes new motion by gyro
disp('Motion caused by gyro causes new motion by gyro');
w = [w1;w2;w3];
wDot = I\(M-cross(w,I*w))