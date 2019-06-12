function [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag, testOmgMag)
%[eigVals rh accels accels2 forces] = StabilityMeasures(.5,.5);
global Xbase
if isempty(Xbase)
    Xbase = findTrim(20);
    Xbase
end

% disp('Baseline');
% [Vcg_dot_base, omg_dot_base, Fsum_base, Msum_base] = testStabilityMeasures(Xbase);

% disp('V_x');
X = Xbase + [0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% % pause
% X = Xout(end,:)';
[Vcg_dot1, omg_dot1, Fsum1, Msum1] = testStabilityMeasures(X);
Vcg_dot1 = Vcg_dot1/testVMag;
omg_dot1 = omg_dot1/testVMag;
% Fsum1 = Fsum1/testVMag;
% Msum1 = Msum1/testVMag;

% Vcg_dot1 = [0 0 0];
% omg_dot1 = [0 0 0];
% Fsum1 = [0 0 0];
% Msum1 = [0 0 0];

% disp('V_y');
X = Xbase + [0; 0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% % pause
% X = Xout(end,:)';
[Vcg_dot2, omg_dot2, Fsum2, Msum2] = testStabilityMeasures(X);
Vcg_dot2 = Vcg_dot2/testVMag;
omg_dot2 = omg_dot2/testVMag;
% Fsum2 = Fsum2/testVMag;
% Msum2 = Msum2/testVMag;

% disp('V_z');
X = Xbase + [0; 0; 0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot3, omg_dot3, Fsum3, Msum3] = testStabilityMeasures(X);
Vcg_dot3 = Vcg_dot3/testVMag;
omg_dot3 = omg_dot3/testVMag;
% Fsum3 = Fsum3/testVMag;
% Msum3 = Msum3/testVMag;

% disp('Omega_x');
X = Xbase + [0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot4, omg_dot4, Fsum4, Msum4] = testStabilityMeasures(X);
Vcg_dot4 = Vcg_dot4/testOmgMag;
omg_dot4 = omg_dot4/testOmgMag;
% Fsum4 = Fsum4/testOmgMag;
% Msum4 = Msum4/testOmgMag;

% disp('Omega_y');
X = Xbase + [0; 0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot5, omg_dot5, Fsum5, Msum5] = testStabilityMeasures(X);
Vcg_dot5 = Vcg_dot5/testOmgMag;
omg_dot5 = omg_dot5/testOmgMag;
% Fsum5 = Fsum5/testOmgMag;
% Msum5 = Msum5/testOmgMag;

% disp('Omega_z');
X = Xbase + [0; 0; 0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot6, omg_dot6, Fsum6, Msum6] = testStabilityMeasures(X);
Vcg_dot6 = Vcg_dot6/testOmgMag;
omg_dot6 = omg_dot6/testOmgMag;
% Fsum6 = Fsum6/testOmgMag;
% Msum6 = Msum6/testOmgMag;

Vcg_dot = [Vcg_dot1; Vcg_dot2; Vcg_dot3; Vcg_dot4; Vcg_dot5; Vcg_dot6];
omg_dot = [omg_dot1; omg_dot2; omg_dot3; omg_dot4; omg_dot5; omg_dot6];
F = [Fsum1; Fsum2; Fsum3; Fsum4; Fsum5; Fsum6];
M = [Msum1; Msum2; Msum3; Msum4; Msum5; Msum6];
% clc
% close all
% disp(['Hover: VDot = ', num2str(Vcg_dot_base), ', OmgDot = ', num2str(omg_dot_base)]);
% disp(['Vx = ', num2str(testVMag), ': VDot = ', num2str(Vcg_dot1), ', OmgDot = ', num2str(omg_dot1)]);
% disp(['Vy = ', num2str(testVMag), ': VDot = ', num2str(Vcg_dot2), ', OmgDot = ', num2str(omg_dot2)]);
% disp(['Vz = ', num2str(testVMag), ': VDot = ', num2str(Vcg_dot3), ', OmgDot = ', num2str(omg_dot3)]);
% disp(['Omgx = ', num2str(testOmgMag), ': VDot = ', num2str(Vcg_dot4), ', OmgDot = ', num2str(omg_dot4)]);
% disp(['Omgy = ', num2str(testOmgMag), ': VDot = ', num2str(Vcg_dot5), ', OmgDot = ', num2str(omg_dot5)]);
% disp(['Omgz = ', num2str(testOmgMag), ': VDot = ', num2str(Vcg_dot6), ', OmgDot = ', num2str(omg_dot6)]);
forces(1:6,1:3) = F;
forces(1:6,4:6) = M;
accels(1:6,1:3) = Vcg_dot;
accels(1:6,4:6) = omg_dot;
forces = forces';
accels = accels';
accels2 = cat(2,accels,[0,-9.81;9.81,0;0,0;0,0;0,0;0,0]);
accels2 = cat(1,accels2,[0,0,0,1,0,0,0,0;0,0,0,0,1,0,0,0]);

accels2(:,6) = [];
accels2(:,3) = [];
accels2(6,:) = [];
accels2(3,:) = [];
eigVals = eig(accels2);
disp(accels2);
disp(eigVals);
rh = [];
% %find equation for whole matrix, fill it in
% % pols = ???
% a = accels2(1,1);
% b = accels2(1,2);
% c = accels2(1,3);
% d = accels2(1,4);
% g = 9.8;
% e = accels2(3,1);
% f = accels2(3,2);
% h = accels2(3,3);
% k = accels2(3,4);
% coef = [ 1, - 2*a - 2*h, a^2 + 4*a*h + b^2 + h^2 + k^2 - 2*c*e + 2*d*f, 2*a*c*e - 2*a*h^2 - 2*a^2*h - 2*b^2*h - 2*a*k^2 - 2*f*g - 2*a*d*f + 2*b*c*f + 2*b*d*e + 2*c*e*h - 2*d*f*h + 2*c*f*k + 2*d*e*k, a^2*h^2 + a^2*k^2 - 2*a*c*e*h - 2*a*c*f*k - 2*a*d*e*k + 2*a*d*f*h + 2*g*a*f + b^2*h^2 + b^2*k^2 + 2*b*c*e*k - 2*b*c*f*h - 2*b*d*e*h - 2*b*d*f*k - 2*g*b*e + c^2*e^2 + c^2*f^2 + d^2*e^2 + d^2*f^2 - 2*g*e*k + 2*g*f*h, 2*b*e*g*h - 2*d*f^2*g - 2*a*f*g*h - 2*d*e^2*g + 2*a*e*g*k + 2*b*f*g*k, e^2*g^2 + f^2*g^2];
% rh = routh_hurwitz(coef);
