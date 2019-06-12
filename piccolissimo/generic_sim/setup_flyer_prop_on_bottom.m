global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m

disp('Setting parameters in setup_flyer_prop_on bottom');
    %Changable Parameters
        %Rotor properties
        h_r = .068; %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .063; %rotor mass in kg
        Ir_r = [203536.211032,87.172143,0.003051;87.172143,26959.742608,-0.127636;0.003051,-0.127636,182477.134918]*1e-9; %Rotor inertia at rotor cg in flyer frame

        %Stator properties
        h_s = -.026; %Height of stator cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_s = .164; %stator mass in kg
        Is_s = [624513.130008,70.172180,379.489186;70.172180,587292.296260,625.057297;379.489186,625.057297,65063.603653]*1e-9; %Stator inertia at stator cg in flyer frame

        %Propeller properties
        h_p = .017; %Height of propeller below ROTOR cg (meters) (positive is down)
        R_p = .178; %Single blade radius (meters)
        beta_p = [0 0 .397 .341 .2998 .258 .224 .194 .171 .151 .136]; %14x6
        chord_p =[.0228 .0171 .00247 .0269 .0299 .0317 .0313 .0291 .0249 .0205 .0146]; %14x6
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)

        %Drag plate properties
        h_d = -.001; %Height of center of drag plate below stator cg (meters) (positive is down)
        R_d1 = .140; %Single plate radius (meters)
        R_d2 = .140; %Single plate radius (meters)
        beta_d = [pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %to be function of r
        chord_d =[.105 .105 .099 .094 .086 .077 .065 .05 .03]; %to be function of r
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)
        
        %Motor properties
        r_m = .26; %ohms
        K_t_m = 1/(740*2*pi/60); %Electromotive force const (K=K_e=K_t=1/K_V) (Nm/A)
        L_m = .2; %Motor inductance (Henry)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.8; %gravity acceleration (m/s^2)
        
        %Simulated Manufacturing Error
        rotorErrorxAngle = deg2rad(0); %2 degrees error
        rotorErroryAngle = 0;
        rotorErrorAngle = [1,0,0;0,cos(rotorErrorxAngle),-sin(rotorErrorxAngle);0,sin(rotorErrorxAngle),cos(rotorErrorxAngle)]*[cos(rotorErroryAngle),0,sin(rotorErroryAngle);0,1,0;-sin(rotorErroryAngle),0,cos(rotorErroryAngle)];