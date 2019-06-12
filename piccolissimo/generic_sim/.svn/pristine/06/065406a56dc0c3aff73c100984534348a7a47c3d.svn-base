global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m

disp('Setting parameters in setup_base_prop_on_bottom');
    %Changable Parameters
        %Rotor properties
        h_r = 0; %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .063; %rotor mass in kg
        Ir_r = [203536.211032,0,0;0,26959.742608,0;0,0,182477.134918]*1e-9; %Rotor inertia at rotor cg in flyer frame

        %Stator properties
        h_s = 0; %Height of stator cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_s = .164; %stator mass in kg
        Is_s = [624513.130008,0,0;0,587292.296260,0;0,0,65063.603653]*1e-9; %Stator inertia at stator cg in flyer frame

        %Propeller properties
        h_p = 0; %Height of propeller below ROTOR cg (meters) (positive is down)
        R_p = .178; %Single blade radius (meters)
        beta_p = 0.2*ones(1,10); %Propeller twist (relative to zero lift)(assumed constant) (radians)
        chord_p = .02*ones(1,10); %Chord lengh (assumed constant) (meters)
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)

        %Drag plate properties
        h_d = 0; %Height of center of drag plate below stator cg (meters) (positive is down)
        R_d1 = .140; %Single plate radius (meters)
        R_d2 = .140; %Single plate radius (meters)
        beta_d = pi/2*ones(1,10); %Propeller twist (relative to zero lift)(assumed constant) (radians)
        chord_d = .15*ones(1,10); %Chord lengh (assumed constant) (meters)
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