
disp('Setting parameters in setup_flyer_test');
clearvars -global
global  Xbase B_p B_d drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m rotorErrorAngle
    Xbase = [3.1426; 0; 0; 0; 0; 0; 0; -297.5032; 37.6379; 0; 0; 0; 0; 0; 0; 0; 0; -2.3980];
    %Changable Parameters
        %Rotor properties
        h_r = 0; %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .055; %rotor mass in kg
        Ir_r = [203536.211032,0,0;0,26959.742608,0;0,0,182477]*1e-9; %Rotor inertia at rotor cg in flyer frame

        %Stator properties
        h_s = 0; %Height of stator cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_s = .156; %stator mass in kg
        Is_s = [624513.130008,0,0;0,587292.296260,0;0,0,65063]*1e-9; %Stator inertia at stator cg in flyer frame

        %Propeller properties
        h_p = 0;%-.018; %-.066 %Height of propeller below ROTOR cg (meters) (positive is down)
        R_p = .178; %Single blade radius (meters)
        beta_p = [0 0 .397 .341 .2998 .258 .224 .194 .171 .151 .136]; %14x6
%         beta_p = interp1(1:length(beta_p),beta_p',1:.1:length(beta_p));
        chord_p =[.0228 .0171 .00247 .0269 .0299 .0317 .0313 .0291 .0249 .0205 .0146]; %14x6
%         chord_p = interp1(1:length(chord_p),chord_p',1:.1:length(chord_p));
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)
        B_p = 2; %number of propeller blades

        %Drag plate properties
        h_d = h_p + 0;%-.010; %Height of center of drag plate below stator cg (meters) (positive is down)
        R_d1 = .150; %Single plate radius (meters)
        R_d2 = .150; %Single plate radius (meters)
        beta_d = [pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %to be function of r
%         beta_d = interp1(1:length(beta_d),beta_d',1:.1:length(beta_d));
%         chord_d =[.105 .105 .099 .094 .086 .077 .065 .05 .03]; %wright bat drag plate
        chord_d = [.05 .05 .05 .05 .05 .05 .05 .05 .05]; %blue drag plate
%         chord_d = interp1(1:length(chord_d),chord_d',1:.1:length(chord_d));
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
        rotorErrorxAngle = deg2rad(1); %2 degrees error
        rotorErroryAngle = deg2rad(1);
        rotorErrorAngle = [1,0,0;0,cos(rotorErrorxAngle),-sin(rotorErrorxAngle);0,sin(rotorErrorxAngle),cos(rotorErrorxAngle)]*[cos(rotorErroryAngle),0,sin(rotorErroryAngle);0,1,0;-sin(rotorErroryAngle),0,cos(rotorErroryAngle)];