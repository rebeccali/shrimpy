disp('Setting parameters in setup_flyer_Spider_V3_centerBatt');
clearvars -global
global R_nu Cd_segmented Cl_segmented Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m rotorErrorAngle Rs_b

    %Changable Parameters
        load('NACA0012_RE360000_segmented_HR.mat');
        
        Xbase = [nan; 0; 0; 0; 0; 0; 0; -553.1437; 30.4464; 0; 0; 0; 0; 0; 0; 0; 0; -5.6941];
        h_cg = -.00244; %cg height from origin
        %% Rotor/Propeller properties
        
        h_r = .02793 - h_cg; %UA hinger rotor %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .01839; %UA hinge rotor %rotor mass in kg
        Ir_r = [ 1715,0,0;0,27649,0;0,0, 27582]*1e-9; %UA hinge rotor %Rotor inertia at rotor cg in flyer frame 
        R_p = .151; %Single blade radius (meters) %UA hinge rotor
        beta_p = atan([0 4.26/23.28 4.31/22.85 3.1/20.88 2.25/18.91 1.55/16.93 .99/14.96 .85/13.48])'; %UA hinge rotor, taken from full blade values at 60-tip every 20mm
        chord_p = [.015 23.67/1000 23.22/1000 21.11/1000 19.04/1000 17/1000 14.99/1000 13.51/1000]'; %UA hinge rotor, taken from full blade values at 60-tip every 20mm

        %% Stator properties
        h_s = -0.0013 - h_cg;
        m_s = .157;%.1183; %stator mass in kg, center battery %.129?
        Is_s = [1110000,0,0;0,1110000,0;0,0,1930000]*1e-9; %Stator inertia at stator cg in flyer frame, center battery
    
        %% Propeller properties
        h_p = .03507 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        beta_p = interp1(1:length(beta_p),beta_p,1:.1:length(beta_p));
        chord_p = interp1(1:length(chord_p),chord_p,1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
        h_d = -.025 - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        R_d1 = .196; %Single plate radius (meters)
        R_d2 = .196; %Single plate radius (meters)
        beta_d = pi-atan([0 0 33.05/23.76 44.39/26.68 44.38/32.7 44.04/41.84 44.27/51.64 44.92/58.79 44.34/63.27 43.17/68.6 43.17/68.6])';
        beta_d = interp1(1:length(beta_d),beta_d,1:.1:length(beta_d));
        chord_d = [0 0 40.71/1000 51.79/1000 55.13/1000 60.75/1000 68.02/1000 73.98/1000 77.26/1000 81.05/1000 81.05/1000]'; 
        chord_d = interp1(1:length(chord_d),chord_d,1:.1:length(chord_d));
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 8; %number of dragplates
        
        %Motor properties
        r_m = .50; %ohms
        K_t_m = 1/(1000*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .02; %Motor inductance (Henry)6
%         r_m = .098; %ohms
%         K_t_m = 1/(2000*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
%         L_m = .01; % intentionally incorrect, integration steps too small to handle dynamics this fast %.000011; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.81; %gravity acceleration (m/s^2)
        
        R_nu = .151;
        
        Xbase(1) = MomentumInflow((m_s+m_r), R_nu, rho);
        
        %Simulated Manufacturing Error
        error_deg = 0;
        rotorErrorxAngle = deg2rad(error_deg); %degrees error
        rotorErroryAngle = deg2rad(0);
        rotorErrorAngle = [1,0,0;0,cos(rotorErrorxAngle),-sin(rotorErrorxAngle);0,sin(rotorErrorxAngle),cos(rotorErrorxAngle)]*[cos(rotorErroryAngle),0,sin(rotorErroryAngle);0,1,0;-sin(rotorErroryAngle),0,cos(rotorErroryAngle)];
        Rs_b = [1,0,0;0,cos(rotorRotationxAngle),-sin(rotorRotationxAngle);0,sin(rotorRotationxAngle),cos(rotorRotationxAngle)]*[cos(rotorRotationyAngle),0,sin(rotorRotationyAngle);0,1,0;-sin(rotorRotationyAngle),0,cos(rotorRotationyAngle)];