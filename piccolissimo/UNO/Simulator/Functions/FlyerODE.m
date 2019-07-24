function dX = FlyerODE(t,X)

global v_clamp r_b nu_style nu Vcg omg omg_r omg_b angles angle_r angle_s pitch1_d F_p M_p F_d M_d i_m v pwm K_t_m r_m Rb_f m g Is_s Ir_r Rr_f I_tot motor_offset
   % (1): nu is the inflow velocity assumed always perpendicular to blade
    % (2): u is the x velocity in the flyer frame
    % (3): v is the y velocity in the flyer frame
    % (4): w is the z velocity in the flyer frame
    % (5): p is the x roll rate in the flyer frame
    % (6): q is the y pitch rate in the flyer frame
    % (7): r is the z yaw rate in the flyer frame
    % (8): omg_r is the rotor yaw rate wrt the stator frame
    % (9): omg_b is the body yaw rate wrt the flyer frame
    % (10): phi is the angle about x from world frame to the flyer frame
    % (11): theta is the angle about y from world frame to the flyer frame
    % (12): psi is the angle about z from world frame to the flyer frame
    % (13): psi_r is the angle about stator frame z to rotor frame
    % (14): psi_s is the angle about flyer frame z to stator frame
    % (15): xe is the world x position
    % (16): ye is the world y position
    % (17): ze is the world z position (+ z is "down")
    % (18): i is the motor current
%         X
        %% State variables
%         nu = X(1);
        Vcg = X(2:4);
        omg = X(5:7);
        omg_r = X(8);
        omg_b = X(9);
        angles = X(10:12);
        angle_r = X(13);
        angle_s = X(14);
    %     xcg = X(15:17);
%         i_m = X(18);

        pitch1_d = 0;

        %% Calculate nu
        if(nu_style == 1)
            %       Numerical fixed point method
            nu = NumericalFixedPointInflow(nu);
            nuDot = 0;
        elseif(nu_style == 0)
            % the simplest way
            nu = X(1);
            nuDot = 0;
        end

        %% Thrust
        [F_p, M_p, F_d, M_d] = ComputeAero(); % Computes aerodynamic forces in the flyer frame
        F_p = sum(F_p,1);
        M_p = sum(M_p,1)+[F_p(3)*motor_offset,0,0]*Rb_f';
        F_d = sum(F_d,1);
        M_d = sum(M_d,1);

        % Rb_f computed in ComputeAero.
        %% Compute motor torque
        %         iDot = (v - K_t_m*(omg_r-omg_b)-r_m*i_m)/L_m;
        iDot = 0;
%         i_m = (v - K_t_m*(omg_r))/r_m; % motor rotational speed is omg_r b/c omg_r is rotor speed wrt stator
        % r_m is motor resistance
        % r_b is battery resistance
        i_m = (v_clamp*pwm-K_t_m*(omg_r))/(pwm*pwm*r_b+r_m);
        v = K_t_m*(omg_r) + i_m*r_m;

        M_mR = [0 0 K_t_m*i_m]; % in rotor frame
        % Rebecca: these two should be the same since this should only be a rotation in the Z axis as long as the prop is aligned with the
        % rotational motion which it always is, so dammit piccoli it's the same thing
        M_mF = M_mR*Rr_f; % in flyer frame

        %% Force caused by gravity in flyer frame
        Fg = [-m*g*sin(angles(2)), m*g*cos(angles(2))*sin(angles(1)), m*g*cos(angles(2))*cos(angles(1))];

        %% CG acceleration in flyer frame
        VcgDot = (1/(m))*(F_p+F_d+Fg)'-cross3(omg,Vcg);

        %% Rotate Inertias
        Ib_bF = Rb_f*Is_s*Rb_f'; %calculated in ComputeAero
        Ir_rF = Rr_f*Ir_r*Rr_f'; %calculated in ComputeAero

        %% Rotor rotational acceleration in rotor frame
        omg_rDot = Ir_r\((M_p*Rr_f'+M_mR).*[0 0 1])';

        %% Body rotational acceleration in flyer frame
        omg_bDot = Ib_bF\((M_d-M_mF).*[0 0 1])'; % z component spins body

        %% CG rotational acceleration
        % (9): omg_b is the body yaw rate wrt the flyer frame
        MGyro1 = -cross3(omg,Ib_bF*(omg+[0; 0; omg_b]));
        MGyro2 = -cross3(omg,Ir_rF*(omg+Rr_f*[0; 0; omg_r+omg_b])); %rotor wrt body, so add ang vels
        omgDot = (Ib_bF+Ir_rF+I_tot)\((M_p+M_d-M_mF.*[1 1 0])'+MGyro1+MGyro2-Ib_bF*omg_bDot-Ir_rF*omg_rDot); %TODO:: omg_rDot in rotor frame, everything else in flyer!!!

        %% World rotational acceleration
        anglesDot = [omg(1)+(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*tan(angles(2)); omg(2)*cos(angles(1))-omg(3)*sin(angles(1));(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*sec(angles(2))];

        %% World velocity
        Rw_f = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        XeDot = Rw_f*Vcg;

        wind_noise_acc = [0 0 0];

   % (1): nu is the inflow velocity assumed always perpendicular to blade % nuDot;
    % (2): u is the x velocity in the flyer frame     % VcgDot(1);
    % (3): v is the y velocity in the flyer frame     % VcgDot(2);
    % (4): w is the z velocity in the flyer frame     % VcgDot(3);
    % (5): p is the x roll rate in the flyer frame     % omgDot(1);
    % (6): q is the y pitch rate in the flyer frame     % omgDot(2);
    % (7): r is the z yaw rate in the flyer frame     % 0;
    % (8): omg_r is the rotor yaw rate wrt the stator frame     % omg_rDot(3);
    % (9): omg_b is the body yaw rate wrt the flyer frame     % omg_bDot(3) + omgDot(3);
    % (10): phi is the angle about x from world frame to the flyer frame     % anglesDot(1);
    % (11): theta is the angle about y from world frame to the flyer frame     % anglesDot(2);
    % (12): psi is the angle about z from world frame to the flyer frame     % anglesDot(3);
    % (13): psi_r is the angle about stator frame z to rotor frame     % omg_r;
    % (14): psi_s is the angle about flyer frame z to stator frame     % omg_b;
    % (15): xe is the world x position     % XeDot(1)+wind_noise_acc(1);
    % (16): ye is the world y position     % XeDot(2)+wind_noise_acc(2);
    % (17): ze is the world z position (+ z is "down")     % XeDot(3)+wind_noise_acc(3);
    % (18): i is the motor current     % iDot];

        dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); 0; omg_rDot(3); omg_bDot(3) + omgDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_b; XeDot(1)+wind_noise_acc(1); XeDot(2)+wind_noise_acc(2); XeDot(3)+wind_noise_acc(3); iDot];

end
