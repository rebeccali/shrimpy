disp('Setting Piccolissimo_V11');
clearvars -global
global motor_offset target_size span_d span_p dSpan_d dSpan_p R_nu Cl_drag Cd_drag Cl_prop Cd_prop Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d m I_tot counter r_m K_t_m L_m Rs_b

% Xbase = [nan; 0; 0; 0; 0; 0; 0; 4129; -269; 0; 0; 0; 0; 0; 0; 0; 0; 0.1851]; %-4.5 deg beta_p %0 deg beta_b %According to findTrim

%0 deg beta_p %0 deg beta_b %According to findTrim
Xbase = [nan;-2.00216238352065e-14;-9.15191250986411e-15;-9.99721101808671e-05;1.85523185398955e-12;-4.85689677629606e-13;0;4940.78051763286;-259.747260994416;-1.98494501464230e-15;2.19425947675265e-15;-8.35391003805421e-27;0;0;2.91828621095746e-15;1.42063992484768e-14;3.17895339412191;0.185100000000000];

%% Heights
h_cg = -2.71/1000; % offset by 3.57 horizontally %cg height from origin (positive is down)
h_r = -3.06/1000 - h_cg; %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
h_d_root = -9.73/1000 - h_cg;
h_d_tip = -2.56/1000 - h_cg;
h_d = h_d_tip; %3/4*h_d_tip+1/4*h_d_root; %Height of 1/4 chord of drag plate below FLYER cg (meters) (positive is down)
h_p = 1.55/1000 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)

% motor_offset = 0;
motor_offset = .008; %.00807;
body_Piccolissimo_V11();
propeller_cheerson_cx10();
motor_cheerson_cx10();
environment_standard;
dSpan = .0005;
slicing_equal;
PiccolissimoControl();