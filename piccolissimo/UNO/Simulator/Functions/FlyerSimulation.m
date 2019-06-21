function [tout, Xout] = FlyerSimulation(X,time_in,argsIn)
    %Syntax:
    %X = [nu; u; v; w; p; q; r; omg_r; omg_b; phi; theta; psi; psi_r; psi_s; xe; ye; ze; i];
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
    %Example:
    %X = [2.912; 0; 0; 0; 0; 0; 0; -278; 45.16; .1; 0; 0; 0; 0; 0; 0; 0; -2.3980]; %14x6
    %or 
    %X = [3.529; 0; 0; 0; 0; 0; 0; -422.6; 38.24; .1; 0; 0; 0; 0; 0; 0; 0; -1.504]; %11x4.7
    %[tout Xout] = yimFlyerLite(X,3,'');
    
    %% Setup
    global pwm pwm_out v_clamp r_b pwm_0 rho chord_d chord_p aoa_p_out aoa_d_out tau_p_out args nu_style counter pitch_p thrust_p_out thrust_d_out RP_tau_p_out RP_tau_d_out i_m_out nu_out time nu h_r m_r Ir_r h_s m_s Is_s span_d span_p S_s_f S_r_f m I_tot r_m K_t_m v v_0 omg_r_0 psi_target psi_des v_out M waypoints waypoint_num target_size MGyros Rs_b flyer_ctrl_fcn
    i_m_out = [];
    thrust_p_out = [];
    thrust_d_out = [];
    RP_tau_p_out = [];
    RP_tau_d_out = [];
    aoa_p_out = [];
    aoa_d_out = [];
    tau_p_out = [];
    pwm_out = [];
    time = time_in(end);
    args = argsIn;

    psi_des = 0;
    v_out = 0;
    M = 0;
    MGyros = 0;

    %Motion Parameters
        waypoint_num = 1;
        psi_target = 0;
        omg_r_0 = X(8);
        v = K_t_m*X(8)+r_m*X(18);
        pwm_sync = roots([X(18)*r_b, -v_clamp, v]);
        pwm_0 = pwm_sync(pwm_sync >=-1 & pwm_sync <=1);
        pwm = pwm_0;
%         v_0 = v;
        nu = X(1);
        if(nnz(strcmpi(args,'simple_nu')))
            nu_style = 0;
        else
            nu_style = 1;
        end
        
    %Calculated Parameters
        S_s_f = [0 0 h_s];
        S_r_f = [0 0 h_r];
        m = m_s+m_r;
        I_tot = (m_s*(S_s_f*(S_s_f')*eye(3)-(S_s_f')*S_s_f)+m_r*(S_r_f*(S_r_f')*eye(3)-(S_r_f')*S_r_f));
        counter = 0;
        nu_out = zeros(length(span_d),1);
        pitch_p = [0 0];
        if(isempty(Rs_b))
            Rs_b = eye(3);
        end
        
        figure(9);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot3(X(15),X(16),X(16),'.','Clipping','off','linewidth',3);
        [a,b,c] = sphere(8);
        axis equal
        a=a*target_size;
        b=b*target_size;
        c=c*target_size;
        for i = 1:size(waypoints,1)
            surf(a+waypoints(i,1),b+waypoints(i,2),c+waypoints(i,3),'FaceColor','none','EdgeAlpha',.2,'Clipping','off');
        end
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        view(3);
        ax = gca;
        ax.CameraUpVector = [0, 0, -1];
%         view(-37.5,-30);
        drawnow();
    tic

    %% Do math
    if length(time_in) == 1
        time_in = linspace(0,time_in,time_in/.001);
    end
    [tout, Xout] = ode113(@FlyerODE,time_in,X,odeset('AbsTol',1e-3,'RelTol',1e-2,'MaxStep',.001,'OutputFcn',flyer_ctrl_fcn));

    %% Make plots
    if (nnz(strcmpi(args,'plot')))
        sz = size(tout);
        orientation = zeros(sz(1),3);
        velocity_world = zeros(sz(1),3);
        for i = 1:sz(1)
            angles(1) = Xout(i,10);
            angles(2) = Xout(i,11);
            angles(3) = Xout(i,12);
            Rw_f = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
            orientation(i,:) = (Rw_f*[0 0 1]')';
            velocity_world(i,:) = Rw_f*Xout(i,2:4)';
        end    
    
        plotFlyer(args, tout, Xout, velocity_world, orientation)
        
    end
end
