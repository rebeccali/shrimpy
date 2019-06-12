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
    % (13): psi_r is the angle about flyer frame z to rotor frame
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
    global time args nu_style counter pitch_p T_p_out T_d_out RP_tau_p_out RP_tau_d_out i_m_out Nu_out nu h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p h_d R_d1 R_d2 beta_d chord_d S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot r_m K_t_m kpv kdv kiv kppsi kdpsi kipsi kpomg kdomg kiomg amplitude V_clamp V V0 omg0 psi0 angleToTarget psidess Vs M waypoints waypoint_num target_size MGyros
    i_m_out = [];
    T_p_out = [];
    T_d_out = [];
    RP_tau_p_out = [];
    RP_tau_d_out = [];
    time = time_in;
    args = argsIn;

    psidess = 0;
    Vs = 0;
    M = 0;
    MGyros = 0;

    %Motion Parameters
%         waypoints = [1 0 0; 1 0 1; 1 1 1; 0 1 1; 0 1 0; 0 0 0];
%         waypoints = [0 0 0; 0 0 1; 1 0 1; 1 0 0; 1 1 0; 1 1 1; 0 1 1; 0 1 0; 0 0 0];
%         waypoints = [10 0 0; 20 0 0];
        waypoints = [0 0 0];
%         waypoints = [999 0 0];

        target_size = .15;
        waypoint_num = 1;
        altitude_des = 0; %meters
        amplitude = .25; %Volt
        angleToTarget = 0;
        V_clamp = 11.1; %Volt
        if nnz(strcmpi(args,'throttle'))
            kpv = 1;%K_t_m;%.1; %1/(1/K_t_m);
            kdv = 0;%K_t_m*.01;
            kiv = 0;%K_t_m;
        else
            kpv = 0;%.1;
            kdv = 0;
            kiv = 0;
        end
        kpomg = 500; %ku = 250 pu = 2.25 seconds when kpv = 1
%         kdomg = 1.2*kpomg/1.6; %10000
        kdomg = -.2*kpomg;
        kiomg = 0;%kpomg*.05;
        kppsi = 0;
        kdpsi = 0;
        kipsi = 0;
        psi0 = -1.3382+0.083393+pi; % phase shift to make vehicle tilt/move in x direction (tilt about -y)
        kpcyc = 1;
        
        omg0 = X(8)-X(9);
        V = K_t_m*(X(8)-X(9))+r_m*X(18);
        V0 = V;
        nu = X(1);
        if(strcmpi(args,'simple_nu'))
            nu_style = 0;
        else
            nu_style = 1;
        end
        
    %Calculated Parameters
        S_s_f = [0 0 h_s];
        S_r_f = [0 0 h_r];
%         area_p = R_p*R_p*pi;
%         area_d = R_d1*R_d1*pi;
        bladeProperties_p = 0;
        bladeProperties_p(1) = h_p;%h_r;
        bladeProperties_p(2) = R_p;
        bladeProperties_p(3) = size(beta_p,2);
        bladeProperties_p = [bladeProperties_p beta_p];
        bladeProperties_p = [bladeProperties_p chord_p];
        bladeProperties_d1 = 0;
        bladeProperties_d1(1) = h_d;%h_s;
        bladeProperties_d1(2) = R_d1;
        bladeProperties_d1(3) = size(beta_d,2);
        bladeProperties_d1 = [bladeProperties_d1 beta_d];
        bladeProperties_d1 = [bladeProperties_d1 chord_d];
        bladeProperties_d2 = 0;
        bladeProperties_d2(1) = h_d;%h_s;
        bladeProperties_d2(2) = R_d2;
        bladeProperties_d2(3) = size(beta_d,2);
        bladeProperties_d2 = [bladeProperties_d2 beta_d];
        bladeProperties_d2 = [bladeProperties_d2 chord_d];
        m = m_s+m_r;
        I_tot = (m_s*(S_s_f*(S_s_f')*eye(3)-(S_s_f')*S_s_f)+m_r*(S_r_f*(S_r_f')*eye(3)-(S_r_f')*S_r_f));
        counter = 0;
        Nu_out = 0;
        pitch_p = [0 0];
        
        figure(9);
        clf;
        hold all;
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
        drawnow();

    tic

    %% Do math
%     odeParams =
%     odeset('AbsTol',1e-8,'RelTol',1e-4,'OutputFcn',@controlLoop);
%     odeParams = odeset('MaxStep',1e-8,'RelTol',1e-4,'OutputFcn',@controlLoop);
%     [tout Xout] = ode113(@yimFlyerLiteODE,linspace(0,time,time/.001),X,odeParams);
%     [tout Xout] = ode45(@yimFlyerLiteODE,linspace(0,time,time/.001),X);
    [tout, Xout] = ode113(@FlyerODE,linspace(0,time_in,time_in/.001),X,odeset('AbsTol',1e-3,'RelTol',1e-2,'MaxStep',.001,'OutputFcn',@FlyerControlLoop));
%     [tout Xout] = ode23tb(@yimFlyerLiteODE,linspace(0,time,time/.001),X,odeset('MaxStep',.001,'OutputFcn',@controlLoop));

    %% Make plots
    figure(1);
    hold off;
    plot(tout,Xout(:,15),'color','r'); %x
    hold on;
    title('Position in World frame');
    plot(tout,Xout(:,16),'color','g'); %y
    plot(tout,Xout(:,17),'color','b'); %z
    
    figure(2); 
    hold off;
    plot(tout,Xout(:,10),'color','r'); %phi
    hold on;
    title('Angle From World to Body');
    plot(tout,Xout(:,11),'color','g'); %theta
    plot(tout,Xout(:,12),'color','b'); %psi
%     plot(tout,Xout(:,13),'color','m'); %rotor psi
%     plot(tout,Xout(:,14),'color','c'); %stator psi
    
    figure(3);
    hold off;
    plot(tout,Xout(:,2),'color','r');%u
    hold on;
    title('Velocity in Flyer Frame');
    plot(tout,Xout(:,3),'color','g');%v
    plot(tout,Xout(:,4),'color','b');%w
    
    figure(4);
    hold off;
    plot(tout,Xout(:,5),'color','r');%p
    hold on;
    title('Angular Velocity in Flyer Frame');
    plot(tout,Xout(:,6),'color','g');%q
    plot(tout,Xout(:,7),'color','b');%r
    plot(tout,Xout(:,8),'color','c');%omg_r
    plot(tout,Xout(:,9),'color','y');%omg_b
    
    figure(5);
    hold off;
%     plot(tout,Xout(:,1)); %nu
    plot(tout(1:length(Nu_out)),Nu_out);
    hold on;
    title('Nu');
    
    sz = size(tout);
    orientation = zeros(sz(1),3);
    for i = 1:sz(1)
        angles(1) = Xout(i,10);
        angles(2) = Xout(i,11);
        angles(3) = Xout(i,12);
        Leb = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        orientation(i,:) = (Leb*[0 0 1]')';
    end
    figure(6);
    hold off;
    plot(tout,orientation(:,1),'color','r');%rotation
    hold on;
    title('Orientation of Spin Axis');
    plot(tout,orientation(:,2),'color','g');%rotation
    plot(tout,orientation(:,3),'color','b');%rotation
    
    figure(7);
    hold off;
%     plot(tout,Xout(:,18),'color','r');
    plot(tout(1:length(i_m_out)),i_m_out,'color','r');
    hold on;
    title('Current and Voltage');
    plot(tout(1:length(Vs)),Vs,'color','g');
    
%     figure(8);
%     hold off;
%     plot(tout(1:end-1),psidess,'color','r');
%     hold on;
%     title('Commanded and error psi');
%     plot(tout,errpsis,'color','g');
%     mean(errpsis)
%     mean(errpsis(10000:end))
%     mean(errpsis(25000:end))
%     
%     figure(9);
%     hold off;
%     plot(tout,vels(:,1),'color','r');
%     hold on;
%     title('World X Y Z velocities');
%     plot(tout,vels(:,2),'color','g');
%     plot(tout,vels(:,3),'color','b');

    figure(9);
    hold off;
%     plot3(Xout(:,15),-Xout(:,16),-Xout(:,17),'Clipping','off','linewidth',3);
    plot3(Xout(:,15),Xout(:,16),Xout(:,17),'Clipping','off','linewidth',3);
%     plot4(Xout(:,15),Xout(:,16),Xout(:,17),tout,'Clipping','off','linewidth',3);
    hold on;
%     title('Position');
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
    hold off;
    z = zeros(length(Xout(:,8)),2);
    plot(tout,Ir_r*(Xout(:,5:7)+[z, Xout(:,8)])',tout,Is_s*(Xout(:,5:7)+[z, Xout(:,9)])',tout,I_tot*Xout(:,5:7)');
    legend('Rx','Ry','Rz','Sx','Sy','Sz','Bx','By','Bz');
    title('Momentum');
    
    figure(12); % Thrust sources
    hold off;
    plot(tout(1:length(T_p_out)),T_p_out, tout(1:length(T_d_out)), T_d_out, tout(1:length(T_d_out)), T_d_out+T_p_out);
    legend('Propeller','Stabilizer','Net');
    title('Thrust');
    
    figure(13); % Thrust sources
    hold off;
    plot(tout(1:length(RP_tau_p_out)),RP_tau_p_out, tout(1:length(RP_tau_d_out)), RP_tau_d_out, tout(1:length(RP_tau_d_out)),RP_tau_p_out+RP_tau_d_out);
    legend('Propeller','Stabilizer','Net');
    title('Torque');
    
    figure(14); % Motor torque
    hold off;
    plot(tout(1:length(i_m_out)),i_m_out*K_t_m);
    title('Motor Torque');