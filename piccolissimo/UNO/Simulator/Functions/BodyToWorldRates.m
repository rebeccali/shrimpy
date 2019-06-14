function [ world_rates ] = BodyToWorldRates( body_rates, world_angles)
%[ world_rates ] = BodyToWorldRates( body_rates, world_angles )
%   Converts local (body) rates to another frame's (world) coordinates by
%   using the angles from the world to body.  Rates are for angles stored
%   as Yaw-Pitch-Roll

world_rates = [body_rates(1)+(body_rates(2)*sin(world_angles(1))+body_rates(3)*cos(world_angles(1)))*tan(world_angles(2)); ...
    body_rates(2)*cos(world_angles(1)) - body_rates(3)*sin(world_angles(1)); ...
    (body_rates(2)*sin(world_angles(1)) + body_rates(3)*cos(world_angles(1)))*sec(world_angles(2))];

end

