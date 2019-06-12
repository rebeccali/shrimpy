disp('Setting parameters in setup_flyer_UNO_v1_centerBatt');
clearvars -global
global  R_nu Cd_segmented Cl_segmented Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m Rs_b

    %Changable Parameters
        load('NACA0012_RE360000_segmented_HR.mat');
        
        
        h_cg = .01112; %cg height from origin
        %% Rotor/Propeller properties
          h_r = .04728 - h_cg; %.0534 11x4.7 %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
          
          %propeller properties 11x4.7
%         Xbase = [nan; 0; 0; 0; 0; 0; 0; -381.7049; 45.7068; 0; 0; 0; 0; 0; 0; 0; 0; -1.242];
%         m_r = .039; %11x4.7 %rotor mass in kg
%         Ir_r = [ 8550,0,0;0,60125,0;0,0, 55551]*1e-9; %11x4.7 no lg %Rotor inertia at rotor cg in flyer frame 
%         R_p = .140; %Single blade radius (meters) %11x4.7
%         beta_p = atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]); %11x4.7
%         chord_p = [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049];%11x4.7
        
        %propeller properties 12x4.5
%         Xbase = [nan; 0; 0; 0; 0; 0; 0; -322.1825; 46.6698; 0; 0; 0; 0; 0; 0; 0; 0; -1.411];
%         m_r = .039; %11x4.7 %rotor mass in kg
%         Ir_r = [ 8550,0,0;0,8300+66500,0;0,0, 8300+66500]*1e-9; %12x4.5 no lg %Rotor inertia at rotor cg in flyer frame 
%         R_p = .150; %Single blade radius (meters) %12x4.5
%         beta_p = atan([0 9.25/13.2 8/22.75 7.4/26.5 5.7/25.5 3.5/20.5 2/11]); %12x4.5
%         chord_p = [.0127 .0161 .0241 .0275 .0261 .0208 .0112]; %12x4.5

        Xbase = [nan; 0; 0; 0; 0; 0; 0; -335.3731; 49.5242; 0; 0; 0; 0; 0; 0; 0; 0; -2.537];
        m_r = .035;
        Ir_r = [ 5179,0,0;0,69000,0;0,0, 71000]*1e-9; % no lg %Rotor inertia at rotor cg in flyer frame 
        R_p = .150; %Single blade radius (meters) %
        beta_p = deg2rad(11)*[1 1 1 1 1 1];
        chord_p = .0195*[1 1 1 1 1 1];

        %% Stator properties
        h_s = .00809 - h_cg;
        m_s = .160; %stator mass in kg
        Is_s = [800000,0,0;0,800000,0;0,0,1280000]*1e-9; %Stator inertia at stator cg in flyer frame % center battery 

    
        %% Propeller properties
        h_p = .04728 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        beta_p = interp1(1:length(beta_p),beta_p',1:.1:length(beta_p));
        chord_p = interp1(1:length(chord_p),chord_p',1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
        h_d = .02 - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        R_d1 = .175; %Single plate radius (meters)
        R_d2 = .175; %Single plate radius (meters)

        span_d = [0 20 51 128.75 175]/1000;
        
        beta_d_nonuniform = deg2rad(180-[90 50 60 27 21]); % UNO V1
        chord_d_nonuniform = [0 .035 .052 .075 .086]; % UNO V1

        beta_d = interp1(span_d,beta_d_nonuniform,linspace(0,max(span_d),100));
        chord_d = interp1(span_d,chord_d_nonuniform,linspace(0,max(span_d),100));

        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 6; %number of dragplates
        
        %Motor properties
        r_m = .50; %ohms
        K_t_m = 1/(1000*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .02; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.8; %gravity acceleration (m/s^2)
        
        R_nu = .16;
        
        Xbase(1) = MomentumInflow((m_s+m_r), R_nu, rho);
        
        %Simulated Manufacturing Error
        rotor_rotation_deg = 0;
        rotorRotationxAngle = deg2rad(rotor_rotation_deg); %degrees error
        rotorRotationyAngle = deg2rad(0);
        Rs_b = [1,0,0;0,cos(rotorRotationxAngle),-sin(rotorRotationxAngle);0,sin(rotorRotationxAngle),cos(rotorRotationxAngle)]*[cos(rotorRotationyAngle),0,sin(rotorRotationyAngle);0,1,0;-sin(rotorRotationyAngle),0,cos(rotorRotationyAngle)];