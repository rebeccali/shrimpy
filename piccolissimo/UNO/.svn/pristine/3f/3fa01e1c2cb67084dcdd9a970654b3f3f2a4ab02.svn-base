disp('Slicing body and blade equally');
%% Nu
R_nu = R_p;
Xbase(1) = MomentumInflow((m_s+m_r), R_nu, rho);

R_max = max(R_d,R_p);
% dSpan = .005;
span = 0:dSpan:R_max;

R_d = R_max;
R_p = R_max;
span_d = span;
span_p = span;
dSpan_d = dSpan;
dSpan_p = dSpan;

beta_p = interp1(span_p_nonuniform,beta_p_nonuniform,span,'linear',0);
chord_p = interp1(span_p_nonuniform,chord_p_nonuniform,span,'linear',0);

beta_d = interp1(span_d_nonuniform,beta_d_nonuniform,span,'linear',0);
chord_d = interp1(span_d_nonuniform,chord_d_nonuniform,span,'linear',0);

chord_p(1:2) = 0; %Why?
chord_d(1:2) = 0; %Why?