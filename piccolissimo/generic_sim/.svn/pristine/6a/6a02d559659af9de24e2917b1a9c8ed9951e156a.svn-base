function [tout, Xout] = yimFlyerLite(X,time,argsIn)
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
    global T_p_out F_p T_d_out F_d RP_tau_p_out M_p RP_tau_d_out M_d i_m i_m_out Cd_segmented Cl_segmented Nu_out Rs_b B_p B_d nu h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot r_m K_t_m L_m kpv kdv kiv kppsi kdpsi kipsi kpomg kdomg kiomg amplitude V_clamp V V0 omg0 psi0 angleToTarget vels errpsis psidess Vs M waypoints waypoint_num target_size MGyros
    
    args = argsIn;
    wind_noise_amplitudes = [0 0 0];
%     vels = zeros(1000*time-1,3);
%     errpsis = zeros(1000*time-1,1); 
    psidess = 0;
    Vs = 0;
    M = 0;
    MGyros = 0;
    %Pick one
%     if isempty(g)
%     setup_base_prop_on_bottom
%         setup_flyer_ideal_prop_on_bottom
%     setup_flyer_prop_on_bottom
%         setup_flyer_ideal_stable
% setup_flyer_V1_5_prop_on_bottom
%         setup_flyerV2_prop_on_bottom
% setup_flyer_test
% setup_flyer_v3
%     end

    %THIS IS A TEMPORARY CHANGE
%     Ir_r = [203536.211032,0,0;0,26959.742608,0;0,0,182477.134918]*1e-9; %Rotor inertia at rotor cg in flyer frame
%     h_d = -h_s;
%     R_d1 = .1665;
%     R_d2 = 0;

%     h_p = -h_p;
%     h_r = -h_r;
%     h_d = -h_d;
%     h_s = -h_s;

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
            kpv = K_t_m;%.1; %1/(1/K_t_m);
            kdv = 0;%K_t_m*.01;
            kiv = 0;%K_t_m;
        else
            kpv = 0;%.1;
            kdv = 0;
            kiv = 0;
        end
        kpomg = 100; %ku = 250 pu = 2.25 seconds when kpv = 1
%         kdomg = 1.2*kpomg/1.6; %10000
        kdomg = 0;%.001*kpomg;
        kiomg = 0;%kpomg*.05;
        kppsi = 0;
        kdpsi = 0;
        kipsi = 0;
        psi0 = -1.3382+0.083393+pi; % phase shift to make vehicle tilt/move in x direction (tilt about -y)
        kpcyc = 1;
        
        omg0 = X(8)-X(9);
%         V = -6.37;
%         V0 = -6.37;
%         V = -4.69;%14x6
%         V0 = -4.69;%14.6
        V = K_t_m*(X(8)-X(9))+r_m*X(18);
        V0 = V;
        nu = X(1);
        
    %Calculated Parameters
        S_s_f = [0 0 h_s];
        S_r_f = [0 0 h_r];
        area_p = R_p*R_p*pi;
        area_d = R_d1*R_d1*pi;
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
%         I_tot = (m_s*(S_s_f')*S_s_f+m_r*(S_r_f')*S_r_f);
        I_tot = (m_s*(S_s_f*(S_s_f')*eye(3)-(S_s_f')*S_s_f)+m_r*(S_r_f*(S_r_f')*eye(3)-(S_r_f')*S_r_f));
        counter = 0;
%         I_tot+Is_s+Ir_r
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
    [tout, Xout] = ode113(@yimFlyerLiteODE,linspace(0,time,time/.001),X,odeset('AbsTol',1e-3,'RelTol',1e-2,'MaxStep',.001,'OutputFcn',@controlLoop));
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
    
    %% Functions
    function wind_noise_acc = GenerateWindNoise(t)
        wind_noise_acc(3) = wind_noise_amplitudes(3)*(sin(t)+cos(17*t))/2;
        wind_noise_acc(1) = wind_noise_amplitudes(1)*(sin(t)+cos(7*t))/2;
        wind_noise_acc(2) = wind_noise_amplitudes(2)*(cos(t)+sin(11*t))/2;
    end

    function stop = controlLoop(t, X, flag) 
        global xOld errvPrev erraltPrev errpsiPrev erraltDotPrev tPrev V_i omg_i psi_i waypoint_old ZPrev reportCounter
        alpha = .2;
        if strcmp(flag,'init')
            V_i = 0;
            omg_i = 0;
            psi_i = 0;
            tPrev = 0;
            erraltPrev = 0;
            erraltDotPrev = 0;
            errvPrev = 0;
            errpsiPrev = 0;
            waypoint_old = 0;
            ZPrev = 0;
            xOld = [];
            reportCounter = 0;
            
        elseif strcmp(flag,'done')

        else
            %% Init normal control function
            stop = false;
            counter = counter + 1;
            reportCounter = reportCounter + 1;
            dt = t(end)-tPrev;
            
            %% Check if z trim is finished
            if(nnz(strcmpi(args,'trim')))
                if ~isempty(xOld)
                    if abs(X(4) - xOld(4)) < .00001 && abs(X(4)) < .01
                        stop = true;
                        omg_r = X(8);
                        omg_b = X(9);
                        i_m = (V - K_t_m*(omg_r-omg_b))/r_m
%                         psidess = psidess(1:counter);
%                         Vs = Vs(1:counter);
                    end
                end
                xOld = X; 
            end
            
            %% Check if we're in a waypoint
            if( (nnz(strcmpi(args,'throttle')) || nnz(strcmpi(args,'cyclic'))) && (sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2 + (waypoints(waypoint_num,3) - X(17))^2) < target_size) )
%                 disp(['At waypoint ' num2str(waypoint_num)]);
                waypoint_num = waypoint_num + 1;
                waypoint_old = 0;
            end
            
            %% Check if we're at the last waypoint
            if (waypoint_num == size(waypoints,1)+1 && nnz(strcmpi(args,'stop')))
                stop = true;
                disp('I made it to my last target');
                return;
            elseif (waypoint_num == size(waypoints,1)+1)
                waypoint_num = 1;
            end

            %% Cyclic control
            amplitude = 0;
            psides = 0;
%             if sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2) > sqrt(target_size^2/2)
                if nnz(strcmpi(args,'cyclic'))
%                     amplitude = min(kpcyc*sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2),30/256*3.7); % proportional pulsing amplitude up to 1
%                     amplitude = .37/2 * min(1,(kpcyc*sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2)) > target_size);% 30/256*4.2
                    amplitude = .37/2;
    %                     amplitude = V_clamp + V0; %because V0 is negative
                    angleToTarget = atan2(waypoints(waypoint_num,2) - X(16),waypoints(waypoint_num,1) - X(15)); % calculate angle to desired waypoint in world
%                     psides = angleToTarget + psi0;
                    
                    % for piccolissimo pulsing
                    angleCurrent = -atan2(X(3),X(2)) - mod(X(12),2*pi); % velocity angle in flyer - world to flyer angle = velocity angle in world
                    angleCurrent = mod(angleCurrent,2*pi); 
                    if angleCurrent > pi
                        angleCurrent = angleCurrent - 2*pi;
                    end
                    angleError = mod(angleToTarget-angleCurrent,2*pi);
                    if angleError > pi
                        angleError = angleError - 2*pi;
                    end
                    % max speed .2 m/s
                    psides = angleToTarget+0.1*angleError + psi0;
                    psides = psi0;
                    
                    psidess(counter) = psides;
                    
                    % Change the blade angle for under actuated simulation
%                     pitch_p = [deg2rad(amplitude)*sin(X(12) + X(13) + psides), -deg2rad(amplitude)*sin(X(12) + X(13) + psides)]; % X(12) = pilot frame psi, X(13) = rotor psi
                end
                
%             end
            
            %% altitude control
            altitude_des = waypoints(waypoint_num,3);
            erralt = altitude_des - X(17,end);
            erraltdot = (1-alpha)*erraltDotPrev + (alpha)*(erraltPrev - X(17,end))/dt;

            omg_i = omg_i + kiomg*erralt*dt;
            omgdes = omg_i + kpomg*erralt + erraltdot*kdomg + omg0;

            errv = omgdes-(X(8)-X(9));
            V_i = V_i + kiv*(errv)*dt;
            V = V_i + kpv*errv + (errv-errvPrev)/dt*kdv + V0; 
            
            %% mix in pulsing and clamp voltage
            V = V + amplitude*sign(cos(X(14)+X(12)+psides));
            if V > V_clamp
                V = V_clamp;
            elseif V < -V_clamp
                V = -V_clamp;
            end

            % Store current values for next calculation
            tPrev = t(end);
            erraltDotPrev = erraltdot;
            erraltPrev = erralt;
            errvPrev = errv;
            waypoint_old = 1;

            % Store any desired variables
            Vs(counter) = V;
            Nu_out(counter) = nu;
            i_m_out(counter) = i_m;
            T_p_out(counter) = F_p(3);
            T_d_out(counter) = F_d(3);
            RP_tau_p_out(counter) = sqrt(M_p(1)*M_p(1) + M_p(2)*M_p(2));
            RP_tau_d_out(counter) = sqrt(M_d(1)*M_d(1) + M_d(2)*M_d(2));

            %Display Simulation Status
%             if mod(counter,1000) == 0
            if toc > 10
%                 atan2(X(3),X(2))
%                 mod(X(12),2*pi)
%                 angleToTarget
%                 angleCurrent
%                 angleError
                figure(9);
            %     plot3(Xout(:,15),-Xout(:,16),-Xout(:,17),'Clipping','off','linewidth',3);
                plot3(X(15),X(16),X(17),'.','Clipping','off','linewidth',3);
            %     plot4(Xout(:,15),Xout(:,16),Xout(:,17),tout,'Clipping','off','linewidth',3);
            %     title('Position');
                drawnow();
                
                disp(['Sim time: ',num2str(t(end),5),', time per calc: ',num2str(toc/reportCounter,5), ', time left: ' num2str((time-t(end))*toc/reportCounter*1000)]);
                tic
                reportCounter = 0;
            end
        end
    end

    function dX = yimFlyerLiteODE(t,X)
%         X
        %% State variables
        nu = X(1);
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
        
        %% Claculate nu
    %   Change in induced air velocity from propeller
        % the right way
%         [~, nuout] = ode45(@inflowODE,linspace(0,1,1/.001),nu,odeset('OutputFcn',@NuControlLoop));
%         nu = nuout(end);
%         nuDot = 0;

%       Numerical fixed point method
        nu = NumericalFixedPointInflow(nu);
        nuDot = 0;
%         
        % the not quite right and potentially unstable way
%         nd1 = -2*nu*abs(-Vcg(3) + nu)/H_p;
%         nd2 = T/(rho*area_p*H_p);
%         nuDot = nd1+nd2;
        
        % the momentum theory way, requires commented nu = X(1);
%         if T < 0
%             nu = 0;
%         else
%             nu = sqrt(T/(2*rho*area_p));
%         end
    
        % the simplest way
%         nu = X(1);
%         nuDot = 0;
        
                %% Thrust
        [F_p, M_p, F_d, M_d] = ComputeAero();
        T = -F_p(3)-F_d(3);
        

        %% Compute motor torque
        %         iDot = (V - K_t_m*(omg_r-omg_b)-r_m*i_m)/L_m;
        
        iDot = 0;
        i_m = (V - K_t_m*(omg_r-omg_b))/r_m;
        
        M_mR = [0 0 K_t_m*i_m]; % in rotor frame
        M_mF = M_mR*Rr_s*Rs_b*Rb_f; % in flyer frame

        %% Force caused by gravity in flyer frame
        Fg = [-m*g*sin(angles(2)), m*g*cos(angles(2))*sin(angles(1)), m*g*cos(angles(2))*cos(angles(1))];

        %% CG acceleration in flyer frame
        VcgDot = (1/(m))*(F_p+F_d+Fg)'-cross3(omg,Vcg);

        %% Rotate Inertias
        Ib_bF = Rb_f*Is_s*Rb_f'; %calculated in ComputeAero
        Ir_rF = Rr_f*Ir_r*Rr_f'; %calculated in ComputeAero

        %% Rotor rotational acceleration in rotor frame
        omg_rDot = Ir_r\((M_p*Rr_f'+M_mR).*[0 0 1])';
%         omg_rDot = Ir_rF\((M_p+M_mR).*[0 0 1])'; %old

        %% Body rotational acceleration in flyer frame
        omg_bDot = Ib_bF\((M_d-M_mF).*[0 0 1])'; % z component spins body

        %% CG rotational acceleration
        MGyro1 = -cross3(omg,Ib_bF*(omg+[0; 0; omg_b]));
        MGyro2 = -cross3(omg,Ir_rF*(omg+Rr_f*[0; 0; omg_r])); 
        omgDot = (Ib_bF+Ir_rF+I_tot)\((M_p+M_d-M_mF.*[1 1 0])'+MGyro1+MGyro2-Ib_bF*omg_bDot-Ir_rF*omg_rDot);

        %% World rotational acceleration
        anglesDot = [omg(1)+(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*tan(angles(2)); omg(2)*cos(angles(1))-omg(3)*sin(angles(1));(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*sec(angles(2))];

        %% World velocity
        Rw_f = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        XeDot = Rw_f*Vcg;

%         wind_noise_acc = GenerateWindNoise(t);
        wind_noise_acc = [0 0 0];

%         dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); omgDot(3); omg_rDot(3); omg_bDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_b; XeDot(1)+wind_noise_acc(1); XeDot(2)+wind_noise_acc(2); XeDot(3)+wind_noise_acc(3); iDot];
        dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); 0; omg_rDot(3); omg_bDot(3) + omgDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_b; XeDot(1)+wind_noise_acc(1); XeDot(2)+wind_noise_acc(2); XeDot(3)+wind_noise_acc(3); iDot];

        
        %% Functions
        function stop = NuControlLoop(~, nu, flag)
            stop = false;
            if strcmp(flag,'init')
            elseif strcmp(flag,'done')
            elseif length(nu)>2
                if abs(nu(end) - nu(end-1)) < .001
                    stop = true;
                end
            end
        end

        function dNu = inflowODE(~,nu_in)
            nu = nu_in;
            [F_p_nu, ~, F_d_nu, ~] = ComputeAero();
            T = -F_p_nu(3)-F_d_nu(3);
            nd1 = -2*nu_in*abs(-Vcg(3) + nu_in)/H_p;
            nd2 = T/(rho*area_d*H_p);
            dNu = nd1+nd2;
        end
        
        function nu_n1 = NumericalFixedPointInflow(nu_n)
            alpha = .1;
            V_inf = sqrt(Vcg(1)*Vcg(1) + Vcg(2)*Vcg(2));
            nu = nu_n;
%             nu_h = sqrt(m*g/(rho*area_p*2));
            [F_p_nu, ~, F_d_nu, ~] = ComputeAero();
            T = -F_p_nu(3)-F_d_nu(3);
            nu_t = sqrt(T/(2*rho*area_d));
            nu_n1 = V_inf + nu_n*(1-alpha) + alpha*(nu_t*nu_t)/sqrt(V_inf*V_inf+nu_n*nu_n);
            
            if(abs((nu_n1-nu_n)/nu_n1) >= .0005)
                nu_n1 = NumericalFixedPointInflow(nu_n1);
            end
        end
        
        function [F_p, M_p, F_d, M_d] = ComputeAero()
            % Computes aerodynamic forces in the flyer frame
            
            % drag blades computed in body frame
            F_d = [0 0 0];
            M_d = [0 0 0];
            for drag_blade = 1:B_d
                % end with 2*pi so future calculations don't care about the blade angle
                blade_angle = 2*pi*drag_blade/B_d; 
                Rb_f = [cos(angle_s + blade_angle), -sin(angle_s + blade_angle), 0; sin(angle_s + blade_angle), cos(angle_s + blade_angle), 0; 0, 0, 1];
                [F, M] = blade(Rb_f*Vcg,Rb_f*(omg+[0; 0; omg_b]),nu,pitch1_d,bladeProperties_d1,rho);
                F_d = F_d + (F*Rb_f);
                M_d = M_d + (M*Rb_f);
            end
            
            % propeller blades computed in rotor frame
            F_p = [0 0 0];
            M_p = [0 0 0];
            for prop_blade = 1:B_p
                % end with 2*pi so future calculations don't care about the blade angle
                blade_angle = 2*pi*prop_blade/B_p;
                Rr_s = [cos(angle_r + blade_angle), -sin(angle_r + blade_angle), 0; sin(angle_r + blade_angle), cos(angle_r + blade_angle), 0; 0, 0, 1];
                Rr_f = Rr_s*Rs_b*Rb_f;
                [F, M] = blade(Rr_f*(Vcg.*[0;0;1]),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),bladeProperties_p,rho); %Vcg.*[0;0;1] for shrowded prop
                F_p = F_p + (F*Rr_f);
                M_p = M_p + (M*Rr_f);
            end
        end

        function [Fb, Mb] = blade(Vcg,omg,nu,pitch1,bladeProperties,rho)
            h = bladeProperties(1);
            R = bladeProperties(2);
            drSteps = bladeProperties(3);
            beta = bladeProperties(4:3+drSteps);
            chord = bladeProperties(4+drSteps:3+drSteps+drSteps);
            dR = R/drSteps;

            Fb = [0 0 0];
            Mb = [0 0 0];
            
            blade_step = 1:drSteps;
            r=dR*(blade_step-1);
            
            %calculate relative wind
            Ut1 = Vcg(1) + omg(2)*h - r*omg(3); %correct
            Up1 = Vcg(3) + r*omg(1) - nu;
            rw = atan2(Up1,Ut1);

            %lift and drag
            aoa1=mod(pi+beta(blade_step)+pitch1+rw,2*pi)-pi;
            ind = int32((aoa1+pi)*100+1);
            qs = rho/2*(Ut1.*Ut1 + Up1.*Up1).*chord(blade_step);
            l1 = qs.*Cl_segmented(ind)';
            d1 = qs.*abs(Cd_segmented(ind))'; 

            %p and t components of force (n is up, t is back)
            n1 = l1.*cos(rw)+d1.*sin(rw);
            t1 = -l1.*sin(rw)+d1.*cos(rw);

            %force vector
            Fb(1) = -sum(t1)*dR;
            Fb(3) = -sum(n1)*dR;

            %moment vector
            Mb(1) = -sum(r.*n1)*dR;
            Mb(2) = -sum(h.*t1)*dR; %correct (sign change for h)
            Mb(3) = sum(r.*t1)*dR;
        end
    end
end

