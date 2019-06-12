testVMag = .4; %m/s
testOmgMag = .5; %rad/s
time = .1;

%setup flyer
% setup_flyer_ideal_prop_on_bottom
% setup_flyer_ideal_stable
setup_flyerV2_prop_on_bottom
h_d = 0;

disp('Baseline');

Xbase = findTrim(10);
% Xbase = [3.020298435969884; 0; 0; 0; 0; 0; 0; -3.197089293498234e+02; 46.753890357773300; 0; 0; 0; 0; 0; 0; 0; 0; -2.787549129133546];
% pause;
[Vcg_dot_base omg_dot_base Fsum_base Msum_base] = testStabilityMeasures(Xbase);

disp('V_x');
X = Xbase + [0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
[Vcg_dot1 omg_dot1 Fsum1 Msum1] = testStabilityMeasures(X);

X = Xbase + [0; testVMag; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
[tout Xout] = yimFlyerLite(X,time,'');
X = Xout(end,:)';
[Vcg_dot2 omg_dot2 Fsum2 Msum2] = testStabilityMeasures(X);

Vcg_dot = [Vcg_dot1; Vcg_dot2];
omg_dot = [omg_dot1; omg_dot2];
F = [Fsum1; Fsum2];
M = [Msum1; Msum2];

ForcesMoments = cell(7);
ForcesMoments(1,2) = java.lang.String('X');
ForcesMoments(1,3) = java.lang.String('Y');
ForcesMoments(1,4) = java.lang.String('Z');
ForcesMoments(1,5) = java.lang.String('L');
ForcesMoments(1,6) = java.lang.String('M');
ForcesMoments(1,7) = java.lang.String('N');
ForcesMoments(2,1) = java.lang.String('before');
ForcesMoments(3,1) = java.lang.String('after');
ForcesMoments(2:3,2:4) = num2cell(F);
ForcesMoments(2:3,5:7) = num2cell(M);


Accelerations = cell(7);
Accelerations(1,2) = java.lang.String('X');
Accelerations(1,3) = java.lang.String('Y');
Accelerations(1,4) = java.lang.String('Z');
Accelerations(1,5) = java.lang.String('L');
Accelerations(1,6) = java.lang.String('M');
Accelerations(1,7) = java.lang.String('N');
Accelerations(2,1) = java.lang.String('before');
Accelerations(3,1) = java.lang.String('after');
Accelerations(2:3,2:4) = num2cell(Vcg_dot);
Accelerations(2:3,5:7) = num2cell(omg_dot);