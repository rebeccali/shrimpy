testVMag = .5; %m/s
testOmgMag = .5; %rad/s
time = .01;

%setup flyer
% setup_flyer_ideal_prop_on_bottom
% setup_flyer_ideal_stable
% setup_flyerV2_prop_on_bottom
% setup_flyer_V1_5_prop_on_bottom

% h_d = -.08892*3/9;

disp('Baseline');
% X = [3.02; 0; 0; 0; 0; 0; 0; -340; 49; 0; 0; 0; 0; 0; 0; 0; 0; -2.524];
% [tout Xout] = yimFlyerLite(X,time);
% Xbase = Xout(end,:)';
if isempty(Xbase)
    Xbase = findTrim(15);
    Xbase
end
% pause;
[Vcg_dot_base omg_dot_base Fsum_base Msum_base] = testStabilityMeasures(Xbase);

disp('V_x');
X = Xbase + [0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% % pause
% X = Xout(end,:)';
[Vcg_dot1 omg_dot1 Fsum1 Msum1] = testStabilityMeasures(X);
% Vcg_dot1 = [0 0 0];
% omg_dot1 = [0 0 0];
% Fsum1 = [0 0 0];
% Msum1 = [0 0 0];

disp('V_y');
X = Xbase + [0; 0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% % pause
% X = Xout(end,:)';
[Vcg_dot2 omg_dot2 Fsum2 Msum2] = testStabilityMeasures(X);

disp('V_z');
X = Xbase + [0; 0; 0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot3 omg_dot3 Fsum3 Msum3] = testStabilityMeasures(X);

disp('Omega_x');
X = Xbase + [0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot4 omg_dot4 Fsum4 Msum4] = testStabilityMeasures(X);

disp('Omega_y');
X = Xbase + [0; 0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot5 omg_dot5 Fsum5 Msum5] = testStabilityMeasures(X);

disp('Omega_z');
X = Xbase + [0; 0; 0; 0; 0; 0; testOmgMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% [tout Xout] = yimFlyerLite(X,time);
% X = Xout(end,:)';
[Vcg_dot6 omg_dot6 Fsum6 Msum6] = testStabilityMeasures(X);

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
ForcesMoments = cell(7);
ForcesMoments(1,2) = java.lang.String('X');
ForcesMoments(1,3) = java.lang.String('Y');
ForcesMoments(1,4) = java.lang.String('Z');
ForcesMoments(1,5) = java.lang.String('L');
ForcesMoments(1,6) = java.lang.String('M');
ForcesMoments(1,7) = java.lang.String('N');
ForcesMoments(2,1) = java.lang.String('u');
ForcesMoments(3,1) = java.lang.String('v');
ForcesMoments(4,1) = java.lang.String('w');
ForcesMoments(5,1) = java.lang.String('p');
ForcesMoments(6,1) = java.lang.String('q');
ForcesMoments(7,1) = java.lang.String('r');
ForcesMoments(2:7,2:4) = num2cell(F);
ForcesMoments(2:7,5:7) = num2cell(M);

Accelerations = cell(7);
Accelerations(1,2) = java.lang.String('X');
Accelerations(1,3) = java.lang.String('Y');
Accelerations(1,4) = java.lang.String('Z');
Accelerations(1,5) = java.lang.String('L');
Accelerations(1,6) = java.lang.String('M');
Accelerations(1,7) = java.lang.String('N');
Accelerations(2,1) = java.lang.String('u');
Accelerations(3,1) = java.lang.String('v');
Accelerations(4,1) = java.lang.String('w');
Accelerations(5,1) = java.lang.String('p');
Accelerations(6,1) = java.lang.String('q');
Accelerations(7,1) = java.lang.String('r');
Accelerations(2:7,2:4) = num2cell(Vcg_dot);
Accelerations(2:7,5:7) = num2cell(omg_dot);

ForcesMoments = ForcesMoments';
Accelerations = Accelerations';
