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
        
        figure(1);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,15),'color','r'); %x
        hold on;
        title('Position in World frame');
        plot(tout,Xout(:,16),'color','g'); %y
        plot(tout,Xout(:,17),'color','b'); %z

        figure(2); 
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,10),'color','r'); %phi
        hold on;
        title('Angle From World to Body');
        plot(tout,Xout(:,11),'color','g'); %theta
        plot(tout,Xout(:,12),'color','b'); %psi
    %     plot(tout,Xout(:,13),'color','m'); %rotor psi
    %     plot(tout,Xout(:,14),'color','c'); %stator psi

        figure(3);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,2),'color','r');%u
        hold on;
        title('Velocity in Flyer Frame');
        plot(tout,Xout(:,3),'color','g');%v
        plot(tout,Xout(:,4),'color','b');%w
        
        figure(33);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,velocity_world(:,1),'color','r');
        hold on;
        title('Velocity in World Frame');
        plot(tout,velocity_world(:,2),'color','g');
        plot(tout,velocity_world(:,3),'color','b');

        figure(4);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,5),'color','r');%p
        hold on;
        title('Angular Velocity in Flyer Frame');
        plot(tout,Xout(:,6),'color','g');%q
        plot(tout,Xout(:,7),'color','b');%r
        plot(tout,Xout(:,8),'color','c');%omg_r
        plot(tout,Xout(:,9),'color','y');%omg_b

        figure(5);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
    %     plot(tout,Xout(:,1)); %nu
        plot(tout(1:size(nu_out,2)),nu_out);
        hold all;
        plot(tout(1:size(nu_out,2)),mean(nu_out,1),'LineWidth',2);
        title('Nu');
        nu_legend = cellstr(num2str(span_d(:)));
        nu_legend{end+1} = 'Mean';
        legend(nu_legend);

        figure(6);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,orientation(:,1),'color','r');%rotation
        hold on;
        title('Orientation of Spin Axis');
        plot(tout,orientation(:,2),'color','g');%rotation
        plot(tout,orientation(:,3),'color','b');%rotation

        figure(7);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
    %     plot(tout,Xout(:,18),'color','r');
        plot(tout(1:length(i_m_out)),i_m_out,'color','r');
        hold on;
        title('Current and Voltage');
        plot(tout(1:length(v_out)),v_out,'color','g');
        plot(tout(1:length(v_out)),v_clamp-pwm_out.*i_m_out*r_b,'color','b');
        legend('i','v_m','v_b');

    %     figure(8);
    %     hold off;
    %     plot(tout(1:end-1),psi_des,'color','r');
    %     hold on;
    %     title('Commanded and error psi');
    %     plot(tout,errpsis,'color','g');
    %     mean(errpsis)
    %     mean(errpsis(10000:end))
    %     mean(errpsis(25000:end))
    %     
        figure(9);
        view(3);
        ax = gca;
        ax.CameraUpVector = [0, 0, -1];
    %     hold off;
    %     plot(tout,vels(:,1),'color','r');
    %     hold on;
    %     title('World X Y Z velocities');
    %     plot(tout,vels(:,2),'color','g');
    %     plot(tout,vels(:,3),'color','b');

    %     figure(9);
    %     hold off;
    % %     plot3(Xout(:,15),-Xout(:,16),-Xout(:,17),'Clipping','off','linewidth',3);
    %     plot3(Xout(:,15),Xout(:,16),Xout(:,17),'Clipping','off','linewidth',3);
    % %     plot4(Xout(:,15),Xout(:,16),Xout(:,17),tout,'Clipping','off','linewidth',3);
    %     hold on;
    % %     title('Position');
    %     [a,b,c] = sphere(8);
    %     axis equal
    %     a=a*target_size;
    %     b=b*target_size;
    %     c=c*target_size;
    %     for i = 1:size(waypoints,1)
    %         surf(a+waypoints(i,1),b+waypoints(i,2),c+waypoints(i,3),'FaceColor','none','EdgeAlpha',.2,'Clipping','off');
    %     end
    %     xlabel('X (m)');
    %     ylabel('Y (m)');
    %     zlabel('Z (m)');

    %     figure(10); % Moment Sources
    %     hold off;
    %     windowSize = 1000;
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,1,1)+M(:,2,1)),'r');
    %     hold on;
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,1,2)+M(:,2,2)),'g');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,3,1)+M(:,4,1)),'b');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,3,2)+M(:,4,2)),'c');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,MGyros(:,1,1)+MGyros(:,2,1)),'y');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,MGyros(:,1,2)+MGyros(:,2,2)),'m');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,1,1)+M(:,2,1)+MGyros(:,1,1)+MGyros(:,2,1)),'k-');
    %     plot(tout(1:end-1),filter(ones(1,windowSize)/windowSize,1,M(:,1,2)+M(:,2,2)+MGyros(:,1,2)+MGyros(:,2,2)),'k:');
    %     legend('Propx','Propy','Dragx','Dragy','Gyrox','Gyroy','x','y');
    %     title('Moments');

        figure(11); % Angular momentum
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        z = zeros(length(Xout(:,8)),2);
        plot(tout,Ir_r*(Xout(:,5:7)+[z, Xout(:,8)])',tout,Is_s*(Xout(:,5:7)+[z, Xout(:,9)])',tout,I_tot*Xout(:,5:7)');
        legend('Rx','Ry','Rz','Sx','Sy','Sz','Bx','By','Bz');
        title('Momentum');

        figure(12); % Thrust sources
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(thrust_p_out)),thrust_p_out, tout(1:length(thrust_d_out)), thrust_d_out, tout(1:length(thrust_d_out)), thrust_d_out+thrust_p_out);
        legend('Propeller','Stabilizer','Net');
        title('Thrust');

        figure(13); % Torque sources
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(RP_tau_p_out)),RP_tau_p_out, tout(1:length(RP_tau_d_out)), RP_tau_d_out, tout(1:length(RP_tau_d_out)),RP_tau_p_out+RP_tau_d_out);
        legend('Propeller','Stabilizer','Net');
        title('Torque');

        figure(14); % Motor torque
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(i_m_out)),i_m_out*K_t_m);
        title('Motor Torque');

        figure(15); % Prop torque
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(i_m_out)),tau_p_out);
        title('Propeller torque directions');

        figure(16);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:size(aoa_p_out,3)),reshape(aoa_p_out(:,1,:),[size(aoa_p_out,1),size(aoa_p_out,3)]));
        title('Propeller 1 aoa vs time');

        figure(17);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:size(aoa_d_out,3)),reshape(aoa_d_out(:,1,:),[size(aoa_d_out,1),size(aoa_d_out,3)]));
        title('Drag 1 aoa vs time');

        figure(18);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_p(chord_p ~= 0),rad2deg(aoa_p_out(chord_p ~= 0,1,end)),'.-')
        title('Propeller 1 aoa vs span (t(end))');
        xlabel('Span (m)');
        ylabel('AOA (deg)');
        grid on;
        ylim([-20 20]);
        xlim([0 max(span_p(end),span_d(end))]);

        figure(19);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_d(chord_d ~= 0),rad2deg(aoa_d_out(chord_d ~= 0,1,end)),'.-');
        title('Stabilizer 1 aoa vs span (t(end))');
        xlabel('Span (m)');
        ylabel('AOA (deg)');
        grid on;
        ylim([-20 20]);
        xlim([0 max(span_p(end),span_d(end))]);

        figure(20);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_d,nu_out(:,end),'.-');
        title('Inflow vs span (t(end))');
        xlabel('Span (m)');
        ylabel('Nu (m/s)');
        grid on;

        figure(21);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        kinematic_viscosity = 1.846e-5/rho; % ~1.568e-5?
        plot(span_d(chord_d ~= 0), abs(Xout(end,9))*span_d(chord_d ~= 0).*chord_d(chord_d ~= 0)/kinematic_viscosity);
        title('Stabilizer Reynolds Number vs Span (t(end))')
        xlabel('Span (m)');
        ylabel('Reynolds Number ()');
        grid on;

        figure(22);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        kinematic_viscosity = 1.846e-5/rho; % ~1.568e-5?
        plot(span_p(chord_p ~= 0), abs(Xout(end,8))*span_p(chord_p ~= 0).*chord_p(chord_p ~= 0)/kinematic_viscosity);
        title('Propeller Reynolds Number vs Span (t(end))')
        xlabel('Span (m)');
        ylabel('Reynolds Number ()');
        grid on;
    end
end

function ClearPlot(do_clear)
    if(do_clear)
        clf;
    else
        hold all;
    end
end