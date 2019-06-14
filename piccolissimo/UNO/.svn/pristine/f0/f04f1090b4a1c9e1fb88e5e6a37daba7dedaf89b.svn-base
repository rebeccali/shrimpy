function [ body_rates ] = WorldToBodyRates( world_rates, world_angles)
%[ body_rates ] = BodyToWorldRates( world_rates, world_angles)
%   Converts world rates to body (local) rates by
%   using the angles from the world to body.  Rates are for angles stored
%   as Yaw-Pitch-Roll

body_rates = [world_rates(1)-world_rates(3)*sin(world_angles(2)); ...
    world_rates(2)*cos(world_angles(1))+world_rates(3)*sin(world_angles(1))*cos(world_angles(2)); ...
    -world_rates(2)*sin(world_angles(1))+world_rates(3)*cos(world_angles(1))*cos(world_angles(2))];

end
