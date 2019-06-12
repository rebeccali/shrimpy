disp('Setting parameters in setup_flyer_WASP_stable');
clearvars -global
global  Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m rotorErrorAngle

    %Changable Parameters
%         Xbase = [4.8883; 0; 0; 0; 0; 0; 0; -549.4082; 34.6068; 0; 0; 0; 0; 0; 0; 0; 0; -2.4730];
%         Xbase = [4.886; 0; 0; 0; 0; 0; 0; -509.7; 42.25; 0; 0; 0; 0; 0; 0; 0; 0; -2.005];
        Xbase = [4.8890; 0; 0; 0; 0; 0; 0; -510.3642; 42.6888; 0; 0; 0; 0; 0; 0; 0; 0; -2.0517];
        h_cg = -.00138; %cg height from origin
        %% Rotor/Propeller properties
        
        h_r = .02623 - h_cg; %.0534 11x4.7 %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .039; %11x4.7 %rotor mass in kg
        Ir_r = [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9; %11x4.7 no lg %Rotor inertia at rotor cg in flyer frame 
        R_p = .140; %Single blade radius (meters) %11x4.7
        beta_p = atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]); %11x4.7
        chord_p = [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049];%11x4.7

        %% Stator properties
        h_s = -0.00464 - h_cg;
        m_s = .32876; %stator mass in kg
        Is_s = [2579947,0,0;0,2593629,0;0,0,4808227]*1e-9; %Stator inertia at stator cg in flyer frame
    
        %% Propeller properties
        h_p = .03995 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        beta_p = interp1(1:length(beta_p),beta_p',1:.1:length(beta_p));
        chord_p = interp1(1:length(chord_p),chord_p',1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
        h_d_from_cg = .00515 %.00515
        h_d = h_d_from_cg - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        R_d1 = .1905; %Single plate radius (meters)
        R_d2 = .1905; %Single plate radius (meters)
%         beta_d = [1 1 1 1 1 1 1 1 1]*deg2rad(180-17); %to be function of r
        beta_d = [1 1 1 1 1 1 1 1 1]*deg2rad(180-60); %to be function of r
        beta_d = interp1(1:length(beta_d),beta_d',1:.1:length(beta_d));
        chord_d = [.00 .00 .06 .06 .06 .06 .06 .06 .06]; %PE wings
        chord_d = interp1(1:length(chord_d),chord_d',1:.1:length(chord_d));
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 8; %number of dragplates
        
        %Motor properties
        r_m = .26; %ohms
        K_t_m = 1/(740*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .02; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.8; %gravity acceleration (m/s^2)
        
        %Simulated Manufacturing Error
        error_deg = 0
        rotorErrorxAngle = deg2rad(error_deg); %degrees error
        rotorErroryAngle = deg2rad(0);
        rotorErrorAngle = [1,0,0;0,cos(rotorErrorxAngle),-sin(rotorErrorxAngle);0,sin(rotorErrorxAngle),cos(rotorErrorxAngle)]*[cos(rotorErroryAngle),0,sin(rotorErroryAngle);0,1,0;-sin(rotorErroryAngle),0,cos(rotorErroryAngle)];