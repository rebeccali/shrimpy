% Test blade math script

Vcg = [0 0 0]';
omg = [0 10 -200]';

setup_flyer_v3

bladeProperties_p = 0;
bladeProperties_p(1) = h_p;
bladeProperties_p(2) = R_p;
bladeProperties_p(3) = size(beta_p,2);
bladeProperties_p = [bladeProperties_p beta_p];
bladeProperties_p = [bladeProperties_p chord_p];

Rr_f1 = rotorErrorAngle*[cos(0), -sin(0), 0; sin(0), cos(0), 0; 0, 0, 1];
[F1, M1] = blade(Rr_f1*Vcg,Rr_f1*(omg),0,0,bladeProperties_p,rho);
F1 = (F1*Rr_f1);
M1 = (M1*Rr_f1);

Rr_f2 = rotorErrorAngle*[cos(0+pi), -sin(0+pi), 0; sin(0+pi), cos(0+pi), 0; 0, 0, 1];
[F2, M2] = blade(Rr_f2*Vcg,Rr_f2*(omg),0,0,bladeProperties_p,rho);
F2 = (F2*Rr_f2);
M2 = (M2*Rr_f2);

F = F1 + F2
M = M1 + M2