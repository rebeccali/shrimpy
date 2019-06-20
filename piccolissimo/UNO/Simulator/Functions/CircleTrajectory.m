function [waypoint] = CircleTrajectory(t)
%% Makes a Cirlce Trajectory for Picolissimo.
%  Inputs:
%  t : a time value of the cirlce trajectory 

posCenter = [998.5, 0, 0];
radius = 0.5;
timeToFlyCircle = 10; % fly circle in 10 seconds

fractionAlongCircle = min(1,max(0,t/timeToFlyCircle));
angleAlongCircle = 2*pi*fractionAlongCircle;

waypoint = [cos(angleAlongCircle), sin(angleAlongCircle), 0]*radius + posCenter;

