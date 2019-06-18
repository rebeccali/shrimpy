function flyerVerticies = getFlyerBodyVerticies()
%% Define Vehicle shape by its verticies
% Define Constants
scaleFactorZ = 1e-3;
scaleFactorY = 17e-3;
scaleFactorX = 4e-3;

bodySFZ = 6e-3;
bodySFY = 2e-3;
bodySFX = 2e-3;

%propeller vertices definiton
x = [0,1,1,0,0,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,0,0,0];
y = [0,0,1,1,0,0;0,1,1,0,0,0;0,1,1,0,1,1;0,0,1,1,1,1];
z = [0,0,0,0,0,1;0,0,0,0,0,1;1,1,1,1,0,1;1,1,1,1,0,1];

%body vertices definiton
bodyX = [0,1,1,0,0,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,0,0,0];
bodyY = [0,0,1,1,0,0;0,1,1,0,0,0;0,1,1,0,1,1;0,0,1,1,1,1];
bodyZ = [0,0,0,0,0,1;0,0,0,0,0,1;1,1,1,1,0,1;1,1,1,1,0,1];

a = 0.5*ones(4,6);
x = x-a;
y = y-a;
z = z-a; 

bodyX = bodyX-a;
bodyY = bodyY-a;
bodyZ = bodyZ-2*a;

z = z.*scaleFactorZ;
y = y.*scaleFactorY;
x = x.*scaleFactorX;

flyerVerticies.bodyZ = bodyZ.*bodySFZ;
flyerVerticies.bodyY = bodyY.*bodySFY;
flyerVerticies.bodyX = bodyX.*bodySFX;

flyerVerticies.propX = x;
flyerVerticies.propY = y;
flyerVerticies.propZ = z;
