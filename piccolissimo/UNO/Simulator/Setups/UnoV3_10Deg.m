disp('Setting UnoV3_10Deg');
clearvars -global
global span_d span_p dSpan_d dSpan_p R_nu Cl_drag Cd_drag Cl_prop Cd_prop Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d m I_tot counter r_m K_t_m L_m Rs_b

%% Heights
h_cg = .03379; %cg height from origin (positive is down)
h_r = 58.09/1000 - h_cg; %.0534 11x4.7 %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
h_d = 0/1000 - h_cg; %Height of 1/4 chord of drag plate below FLYER cg (meters) (positive is down)
h_p = 67.22/1000 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)

body_UnoV3_3RimBatteries();
propeller_LLFlap_10_deg_0_twist_Quanum2206();
motor_MT2206_2000KV();
environment_standard;
slicing_equal;
UnoControl();