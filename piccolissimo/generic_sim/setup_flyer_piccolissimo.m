disp('Setting parameters in setup_flyer_piccolissimo.m');
clearvars -global
global Cd_segmented Cl_segmented Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m Rs_b

    %Changable Parameters
        load('Cd_segmented_lowRe');
        load('Cl_segmented_lowRe');

        Xbase = [3.92; 0; 0; 0; 0; 0; 0; 5420.1; -399.79; 0; 0; 0; 0; 0; 0; 0; 0; 0.4707]; %8 deg ang %According to findTrim

        h_cg = -.00347; %cg height from origin
        %% Rotor/Propeller properties
        
        h_r = .0027 - h_cg; % Cheerson %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .00007; % Cheerson %rotor mass in kg
        Ir_r = [ .1,0,0;0,2.71,0;0,0, 2.73]*1e-9; %Cheerson %Rotor inertia at rotor cg in flyer frame 
        R_p = .0146; %Cheerson %Single blade radius (meters) 
        beta_p = -atan([0 .56/3.17 .72/3.91 .75/4.06 .63/3.48 .15/.89])'; %Cheerson taken every 2.92mm
        chord_p = [0 3.22/1000 3.98/1000 4.12/1000 3.54/1000 .90/1000]'; %Cheerson taken every 2.92mm

        %% Stator properties
        h_s = -0.00357 - h_cg;
        m_s = .00486-m_r;%.00454; %stator mass in kg
        Is_s = [591,0,0;0,600,0;0,0,960]*1e-9; %Stator inertia at stator cg in flyer frame
    
        %% Propeller properties
        h_p = .0027 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        beta_p = interp1(1:length(beta_p),beta_p,1:.1:length(beta_p));
        chord_p = interp1(1:length(chord_p),chord_p,1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
        h_d = -.007 - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        R_d1 = .02; %Single plate radius (meters)
        R_d2 = .02; %Single plate radius (meters)
        beta_d = deg2rad([0 81 72 63 54 48 43 38])'; %
        beta_d = interp1(1:length(beta_d),beta_d,1:.1:length(beta_d));
        chord_d = [0 13.5 12.825 12.15 11.475 9.45 7.425 4.05]/1000'; 
        chord_d = interp1(1:length(chord_d),chord_d,1:.1:length(chord_d));
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 6; %number of dragplates
        
        %Motor properties
        r_m = 1.97; %ohms
        K_t_m = 1/(20000*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .01; %cheerson 6uf % intentionally incorrect, integration steps too small to handle dynamics this fast %.000011; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.81; %gravity acceleration (m/s^2)
        
        %Simulated Manufacturing Error
        rotor_rotation_deg = 8
        rotorRotationxAngle = deg2rad(rotor_rotation_deg); %degrees error
        rotorRotationyAngle = deg2rad(0);
        Rs_b = [1,0,0;0,cos(rotorRotationxAngle),-sin(rotorRotationxAngle);0,sin(rotorRotationxAngle),cos(rotorRotationxAngle)]*[cos(rotorRotationyAngle),0,sin(rotorRotationyAngle);0,1,0;-sin(rotorRotationyAngle),0,cos(rotorRotationyAngle)];