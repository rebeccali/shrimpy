disp('Setting parameters in setup_flyer_V3');
clearvars -global
global  Xbase B_p B_d h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m rotorErrorAngle
    %Changable Parameters
    Xbase = [3.1426; 0; 0; 0; 0; 0; 0; -297.5032; 37.6379; 0; 0; 0; 0; 0; 0; 0; 0; -2.3980];
        %% Rotor/Propeller properties
%         h_cg = -.01553; %14x6 from origin
%         h_r = .03329 - h_cg;%14x6 
%         m_r = .05527; %14x6 
%         Ir_r = [ 14895,0,0;0,191513,0;0,0,181538 ]*1e-9; %14x6 no lg
%         R_p = .178; %Single blade radius (meters) %14x6
%         beta_p = [0 0 .397 .341 .2998 .258 .224 .194 .171 .151 .136]; %14x6
%         chord_p =[.0228 .0171 .00247 .0269 .0299 .0317 .0313 .0291 .0249 .0205 .0146]; %14x6

        h_cg = -.02074; %11x4.7 from origin
        h_r = .02739 - h_cg; %.0534 11x4.7 %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .039; %11x4.7 %rotor mass in kg
        Ir_r = [ 10897.80,0,0;0,62472.25,0;0,0,55551.13]*1e-9; %11x4.7 no lg %Rotor inertia at rotor cg in flyer frame 
        R_p = .140; %Single blade radius (meters) %11x4.7
        beta_p = atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]); %11x4.7
        chord_p = [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049];%11x4.7

        %% Stator properties
        h_s = -.03278 - h_cg;
        m_s = .156; %stator mass in kg
        Is_s = [ 234954.16,0,0;0,162086.24,0;0,0,107002.91]*1e-9; %Stator inertia at stator cg in flyer frame
    
        %% Propeller properties
        h_p = .043 - h_cg; %Height of propeller below FLYER cg (meters) (positive is down)
        beta_p = interp1(1:length(beta_p),beta_p',1:.1:length(beta_p));
        chord_p = interp1(1:length(chord_p),chord_p',1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades
        
        %Drag plate properties
%         h_d = -.024 - h_cg; %Height of center of drag plate below FLYER cg (meters) (positive is down) default: -.014 - h_cg;
        h_d_from_cg = -.054
        h_d = h_d_from_cg - h_cg; %Test configs, COP=COM at ~-.029, with 5 degrees -.049 is best

        R_d1 = .150; %Single plate radius (meters)
        R_d2 = .150; %Single plate radius (meters)
        beta_d = [pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %to be function of r
        beta_d = interp1(1:length(beta_d),beta_d',1:.1:length(beta_d));
%         chord_d =[.105 .105 .099 .094 .086 .077 .065 .05 .03]; %wright bat drag plate
        chord_d = [.05 .05 .05 .05 .05 .05 .05 .05 .05]; %blue drag plate
        chord_d = interp1(1:length(chord_d),chord_d',1:.1:length(chord_d));
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        B_d = 2; %number of dragplates
        
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