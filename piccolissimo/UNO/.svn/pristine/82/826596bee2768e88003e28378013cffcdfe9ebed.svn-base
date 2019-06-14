disp('Setting parameters in setup_flyer_piccolissimo_V7.m');
clearvars -global
global motor_offset dSpan_p dSpan_d Cd_prop Cl_prop Cd_drag Cl_drag span_p span_d R_nu Cd_segmented Cl_segmented Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m Rs_b

    %Changable Parameters
        Cld = load('FlatPlateLentink_RE14000_segmented_HR.mat');
        Cl_drag = Cld.Cl_segmented;
        Cd_drag = Cld.Cd_segmented;
        Cl_prop = Cld.Cl_segmented;
        Cd_prop = Cld.Cd_segmented;

%         Xbase = [nan; 0; 0; 0; 0; 0; 0; 5198; -306.2; 0; 0; 0; 0; 0; 0; 0; 0; 0.3727]; %According to findTrim
%         Xbase = [nan; 0; 0; 0; 0; 0; 0; 4897; -303.1; 0; 0; 0; 0; 0; 0; 0; 0; 0.362]; %-1 deg beta_p %According to findTrim
%         Xbase = [nan; 0; 0; 0; 0; 0; 0; 4206; -300.2; 0; 0; 0; 0; 0; 0; 0; 0; 0.3525]; %-4 deg beta_p %According to findTrim
        Xbase = [nan; 0; 0; 0; 0; 0; 0; 4119; -296.8; 0; 0; 0; 0; 0; 0; 0; 0; 0.3541]; %-4.5 deg beta_p %.5 deg beta_b %According to findTrim

        h_cg = -.00281; %cg height from origin
        %% Rotor/Propeller properties
        
        h_r = .0027 - h_cg; % Cheerson %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .00007 + .00025; % Cheerson %rotor mass in kg
        Ir_r = [ .1+1.36,0,0;0,2.71+1.36,0;0,0, 2.73+0.97]*1e-9; %Cheerson prop plus motor rotor %Rotor inertia at rotor cg in flyer frame 
        R_p = .0146; %Cheerson %Single blade radius (meters) 
        span_p_nonuniform = 0:.00292:R_p;
        beta_p_nonuniform = -deg2rad(4.5)-atan([0 .56/3.17 .72/3.91 .75/4.06 .63/3.48 .15/.89])'; %Cheerson taken every 2.92mm
        chord_p_nonuniform = [0 3.22/1000 3.98/1000 4.12/1000 3.54/1000 .90/1000]'; %Cheerson taken every 2.92mm

        %% Stator properties
        h_s = -0.00291 - h_cg;
        m_s = .00413-m_r;%.00454; %stator mass in kg
        Is_s = [463,0,0;0,468,0;0,0,783]*1e-9; %Stator inertia at stator cg in flyer frame
    
        %% Propeller properties
        h_p = .004 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
        h_d = -.001 - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        R_d = .02; %Single plate radius (meters)
        span_d_nonuniform = 0:.004:R_d;
        beta_d_nonuniform = deg2rad(.5)+atan([0 2.51/.68 5.3/2.18 9.52/6.08 6.83/5.58 4.27/3.85])'; % Taken every 4mm 
        chord_d_nonuniform = [0 2.6 5.74 11.36 8.88 5.81]/1000'; % Taken every 4mm
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 6; %number of dragplates
        
        %Motor properties
        r_m = 1.97; %ohms
        K_t_m = 1/(20500*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .01; %cheerson 6uf % intentionally incorrect, integration steps too small to handle dynamics this fast %.000011; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.81; %gravity acceleration (m/s^2)
        
        dSpan = .0005;
        slicing_equal;
        
        PiccolissimoControl;
        
        %Motor Offset
        motor_offset = .004;
%         motor_offset = 0;
        
        %Motor Angle
        rotor_rotation_deg = 0;
        rotorRotationxAngle = deg2rad(rotor_rotation_deg); %degrees error
        rotorRotationyAngle = deg2rad(0);
        Rs_b = [1,0,0;0,cos(rotorRotationxAngle),-sin(rotorRotationxAngle);0,sin(rotorRotationxAngle),cos(rotorRotationxAngle)]*[cos(rotorRotationyAngle),0,sin(rotorRotationyAngle);0,1,0;-sin(rotorRotationyAngle),0,cos(rotorRotationyAngle)];