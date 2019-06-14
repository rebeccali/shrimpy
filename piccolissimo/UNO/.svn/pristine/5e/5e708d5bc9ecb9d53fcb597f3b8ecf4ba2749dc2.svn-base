disp('Setting UnoOffset_actual_gemfan_10x4_5');
clearvars -global
global motor_offset target_size span_d span_p dSpan_d dSpan_p R_nu Cl_drag Cd_drag Cl_prop Cd_prop Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d m I_tot counter r_m K_t_m L_m Rs_b

%% Heights
h_cg = .01999; %cg height from origin (positive is down)
h_r = 53.49/1000 - h_cg; %.0534 11x4.7 %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
h_d = 0/1000 - h_cg; %Height of 1/4 chord of drag plate below FLYER cg (meters) (positive is down)
h_p = 60.71/1000 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)

motor_offset = 0;
body_UnoOffset_actual();
propeller_gemfan_10x4_5();
motor_MT2206_2000KV();
environment_standard;
dSpan = .005;
slicing_equal;
UnoControl();