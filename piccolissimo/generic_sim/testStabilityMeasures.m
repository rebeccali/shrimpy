function [Vcg_dot omg_dot meansF meansM] = testStabilityMeasures(X)

global V rotorErrorAngle B_p B_d drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m
    persistent Cd_segmented Cl_segmented
    if isempty(Cd_segmented)
        load('Cd_segmented');
        load('Cl_segmented');
    end
%setup velocities, etc.
    angle_r_orig = 0;%X(13);
    angle_s_orig = 0;%X(14);
    S_s_f = [0 0 h_s];
    S_r_f = [0 0 h_r];
    area_p = R_p*R_p*pi;
%         area_d = R_d*R_d*pi;
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
    
    
    %setup blade orientations
    rotation_positions = 24;
    rotation_step = 2*pi/rotation_positions; %only do half of a rotation
    k = 1;
    for i = 0:(rotation_positions-1)
        for j = 0:(rotation_positions-1)
            X(13) = angle_r_orig + i*rotation_step;
            X(14) = angle_s_orig + j*rotation_step;
            %run test
            [dX(k,:) Fsum(k,:) Msum(k,:)] = FlyerODE(0,X);
            k = k+1;
        end
    end
    %sum tests
    means = mean(dX,1);
    meansF = mean(Fsum,1);
    meansM = mean(Msum,1);
    Vcg_dot = means(2:4);
    omg_dot = means(5:7);

    %% Functions
    function stop = controlLoop(t, X, flag) 
        global xOld errvPrev erraltPrev errpsiPrev tPrev V_i omg_i psi_i waypoint_old ZPrev MGyro1 MGyro2
        if strcmp(flag,'init')
            V_i = 0;
            omg_i = 0;
            psi_i = 0;
            tPrev = 0;
            erraltPrev = 0;
            errvPrev = 0;
            errpsiPrev = 0;
            waypoint_old = 0;
            ZPrev = 0;
            xOld = [];
        elseif strcmp(flag,'done')

        else
            stop = false;
            counter = counter + 1;
            if(nnz(strcmpi(args,'trim')))
                if ~isempty(xOld)
                    if abs(X(4) - xOld(4)) < .000001 && abs(X(4)) < .001
                        stop = true;
                        psidess = psidess(1:counter);
                        Vs = Vs(1:counter);
                        M = M(1:counter,:,:);
                        MGyros = MGyros(1:counter,:,:);
                    end
                end
                xOld = X; 
            end
            if((nnz(strcmpi(args,'throttle')) || nnz(strcmpi(args,'cyclic')) && sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2 + (waypoints(waypoint_num,3) - X(17))^2) < target_size))
                if(t(end) > 3)
    %                 disp('I am at target: ');
    %                 waypoint_num
                    waypoint_num = waypoint_num + 1;
                    waypoint_old = 0;
                end
            end
            if (waypoint_num == size(waypoints,1)+1 && nnz(strcmpi(args,'stop')))
                stop = true;
                disp('I made it to my last target');
            elseif (waypoint_num == size(waypoints,1)+1)
                waypoint_num = 1;
            end
            dt = t(end)-tPrev;

            altitude = waypoints(waypoint_num,3);
            amplitude = 0;
            if sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2) > sqrt(target_size^2/2)
                if nnz(strcmpi(args,'cyclic'))
                    amplitude = min(sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2),3);
    %                     amplitude = V_clamp + V0; %because V0 is negative
                end
                translation_angle = atan2(waypoints(waypoint_num,2) - X(16),waypoints(waypoint_num,1) - X(15));
            end
            psides = translation_angle + psi0;
            psidess(counter) = psides;

            erralt = altitude - X(17,end);
            erraltdot = (erraltPrev - X(17,end))/dt;

            omg_i = omg_i + kiomg*erralt*dt;
            omgdes = omg_i + kpomg*erralt + erraltdot*kdomg + omg0;

            errv = omgdes-(X(8)-X(9));
            V_i = V_i + kiv*(errv)*dt;
            V = V_i + kpv*errv + (errv-errvPrev)/dt*kdv + V0;
    %             if sin(X(14)-X(12)+psides) > 0
    %                 V = V + amplitude;
    %             else
    %                 V = V - amplitude;
    %             end
    %             V = V + amplitude*sin(X(14)-X(12)+psides);
            V = V - amplitude*cos(X(14)-X(12)+psides);
            if V > V_clamp
                V = V_clamp;
            elseif V < -V_clamp
                V = -V_clamp;
            end

            erraltPrev = erralt;
            errvPrev = errv;
            waypoint_old = 1;

            % Store any desired variables
            Vs(counter) = V;
    %         M(counter,1,:) = M1;
    %         M(counter,2,:) = M2;
    %         M(counter,3,:) = M3;
    %         M(counter,4,:) = M4;
%             MGyros(counter,1,:) = MGyro1;
%             MGyros(counter,2,:) = MGyro2;

            %Display Simulation Status
            if mod(counter,1000) == 0
                disp(['Sim time: ',num2str(t(end),5),', time per calc: ',num2str(toc/1000,5)]);
    %             moms = M1+M2+M3+M4;
    %             
    %             
    %             omg = X(5:7);
    %             omg_r = X(8);
    %             omg_s = X(9);
    %             % Need to low pass this for it to make any sense!!!
    %             gyromoms = -cross3(omg,Is_sF*(omg+[0; 0; omg_s]))'-cross3(omg,Ir_rF*(omg+[0; 0; omg_r]))';
    %             disp(['aero mom: ', num2str(moms)]);
    %             disp(['gyro mom: ', num2str(gyromoms)]);
    %             disp(['Net aero xy mom: ', num2str(sqrt(moms(1)*moms(1) + moms(2)*moms(2)))]);
    %             disp(['Net gyro xy mom: ', num2str(sqrt(gyromoms(1)^2 + gyromoms(2)^2)) char(10)]);
                tic
            end
        end
    end

    function [dX, Fsum, Msum] = FlyerODE(t,X)

    % global  Ir_r Is_s H_p rho g area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot V r_m K_t_m L_m nu rotorErrorAngle Is_sF Ir_rF MGyro1 MGyro2

        %State variables
        nu = X(1);
        Vcg = X(2:4);
        omg = X(5:7);
        omg_r = X(8);
        omg_s = X(9);
        angles = X(10:12);
        angle_r = X(13);
        angle_s = X(14);
    %     xcg = X(15:17);
        i_m = X(18);

        pitch1_p = 0;
        pitch2_p = 0;
        pitch1_d = 0;
        pitch2_d = 0;

        %Thrust
        [F_p, M_p, F_d, M_d] = ComputeAero();
        Fsum = F_p + F_d;
        Msum = M_p + M_d;   
%         disp(F_d)

        T = -F_p(3)-F_d(3);
    %   Change in induced air velocity from propeller
    % TODO::I could use a better algorithm
%         nd1 = -2*nu*abs(-Vcg(3) + nu)/H_p;
%         nd2 = T/(rho*area_p*H_p);
%         nuDot = nd1+nd2;
    %     nu = sqrt(T/(2*rho*area_p));
        nuDot = 0;

        %Compute motor torque
        M_m = [0 0 K_t_m*i_m];
        iDot = (V - K_t_m*(omg_r-omg_s)-r_m*i_m)/L_m;

        %Force caused by gravity
        Fg = [-m*g*sin(angles(2)), m*g*cos(angles(2))*sin(angles(1)), m*g*cos(angles(2))*cos(angles(1))];

        %CG acceleration
        VcgDot = (1/(m))*(F_p+F_d+Fg)'-cross3(omg,Vcg);

        %Rotate Inertias
        Is_sF = Rs_f*Is_s*Rs_f'; %calculated in ComputeAero
        Ir_rF = Rr_f*Ir_r*Rr_f'; %calculated in ComputeAero

        %Rotor rotational acceleration
        omg_rDot = Ir_rF\((M_p+M_m).*[0 0 1])';

        %Stator rotational acceleration
        omg_sDot = Is_sF\((M_d-M_m).*[0 0 1])';

        %CG rotational acceleration
        MGyro1 = -cross3(omg,Is_sF*(omg+[0; 0; omg_s]));
        MGyro2 = -cross3(omg,Ir_rF*(omg+[0; 0; omg_r]));
        omgDot = (Is_sF+Ir_rF+I_tot)\((M_p+M_d)'+MGyro1+MGyro2-Is_sF*omg_sDot-Ir_rF*omg_rDot);

        %World rotational acceleration
        anglesDot = [omg(1)+(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*tan(angles(2)); omg(2)*cos(angles(1))-omg(3)*sin(angles(1));(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*sec(angles(2))];

        %World velocity
        Ri_f = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        XeDot = Ri_f*Vcg;

        dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); omgDot(3); omg_rDot(3); omg_sDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_s; XeDot(1); XeDot(2); XeDot(3); iDot];

        function [F_p, M_p, F_d, M_d] = ComputeAero()
            F_p = [0 0 0];
            M_p = [0 0 0];
            for prop_blade = 1:B_p
                blade_angle = 2*pi*prop_blade/B_p;
%                 rotor_error_angle = rotorErrorAngle*[cos(angle_s), -sin(angle_s), 0; sin(angle_s), cos(angle_s), 0; 0, 0, 1]; %body spins so the error spins
                Rr_f = [cos(angle_r + blade_angle), -sin(angle_r + blade_angle), 0; sin(angle_r + blade_angle), cos(angle_r + blade_angle), 0; 0, 0, 1];
                [F, M] = blade(Rr_f*Vcg,Rr_f*(omg+[0; 0; omg_r]),nu,pitch1_p,bladeProperties_p,rho);
                F_p = F_p + (F*Rr_f);
                M_p = M_p + (M*Rr_f);
            end

            F_d = [0 0 0];
            M_d = [0 0 0];
            for drag_blade = 1:B_d
                blade_angle = 2*pi*drag_blade/B_d;
                Rs_f = [cos(angle_s + blade_angle), -sin(angle_s + blade_angle), 0; sin(angle_s + blade_angle), cos(angle_s + blade_angle), 0; 0, 0, 1];
                [F, M] = blade(Rs_f*Vcg,Rs_f*(omg+[0; 0; omg_s]),nu,pitch1_d,bladeProperties_d1,rho);
                F_d = F_d + (F*Rs_f);
                M_d = M_d + (M*Rs_f);
            end
        end

        function [Fb, Mb] = blade(Vcg,omg,nu,pitch1,bladeProperties,rho)
%             Vcg
            h = bladeProperties(1);
            R = bladeProperties(2);
            drSteps = bladeProperties(3);
            beta = bladeProperties(4:3+drSteps);
            chord = bladeProperties(4+drSteps:3+drSteps+drSteps);

            Fb = [0 0 0];
            Mb = [0 0 0];

            for blade_step = 1:drSteps
                r=R*(blade_step-1)/drSteps;
                %calculate relative wind
            %     Ut1 = Vcg(1) - (omg(2)*h+r*omg(3)); %old
                Ut1 = Vcg(1) + omg(2)*h - r*omg(3); %correct
                Up1 = Vcg(3) + r*omg(1) - nu; % Up1 is -z

                %lift and drag
                aoa1=mod(pi+beta(blade_step)+pitch1+atan2(Up1,Ut1),2*pi)-pi;
                l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(blade_step)*Cl_segmented(int32((aoa1+pi)*100+1));
                d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(blade_step)*abs(Cd_segmented(int32((aoa1+pi)*100+1)));
%                 l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(blade_step)*abs(Cl_segmented(int32((aoa1+pi)*100+1)));
%                 d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(blade_step)*abs(Cd_segmented(int32((aoa1+pi)*100+1)));
            %     l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa1));
            %     d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa1));

                %p and t components of force (n is up, t is back)
                n1 = l1*cos(atan2(Up1,Ut1))+d1*sin(atan2(Up1,Ut1));
                t1 = -l1*sin(atan2(Up1,Ut1))+d1*cos(atan2(Up1,Ut1));

                %force vector
                Fb(1) = Fb(1)-t1*(R/drSteps); %- force felt backwards to orient forwards
                Fb(3) = Fb(3)-n1*(R/drSteps); %- force felt up to orient down

                %moment vector
                Mb(1) = Mb(1) - r*n1*(R/drSteps); 
                Mb(2) = Mb(2) - h*t1*(R/drSteps); %correct (sign change for h)
                Mb(3) = Mb(3) + r*t1*(R/drSteps);
            end
        end
    end

end

% function [dX Fsum Msum] = FlyerODE(t,X)
%     
% global  V drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter r_m K_t_m L_m
%     
%     %State variables
%     nu = X(1);
%     Vcg = X(2:4);
%     omg = X(5:7);
%     omg_r = X(8);
%     omg_s = X(9);
%     angles = X(10:12);
%     angle_r = X(13);
%     angle_s = X(14);
%     xcg = X(15:17);
%     i_m = X(18);
% 
%     pitch1_p = 0;
%     pitch2_p = 0;
%     pitch1_d = 0;
%     pitch2_d = 0;
%     
%     [tout nuout] = ode113(@inflowODE,linspace(0,10,10/.001),nu,odeset('MaxStep',.001,'OutputFcn',@controlLoop));
%     nu = nuout(end);
%     function stop = controlLoop(t, nu, flag)
%         global lastNu;
%         stop = false;
%         if strcmp(flag,'init')
%             lastNu = nu;
%         elseif strcmp(flag,'done')
%             
%         else
%             if abs(nu - lastNu) < .0001
%                 stop = true;
%             end
%             lastNu = nu;
%         end
%     end
%     
%     function dNu = inflowODE(t,nu)
%         %TODO::use me in yimFlyerLiteODE down there \/
%         %Thrust
%         Rr_f1 = [cos(angle_r), -sin(angle_r), 0; sin(angle_r), cos(angle_r), 0; 0, 0, 1];
%         [F1 M1] = blade(Rr_f1*Vcg,Rr_f1*(omg+[0; 0; omg_r]),nu,pitch1_p,bladeProperties_p,rho);
%         F1 = (F1*Rr_f1);
% 
%         Rr_f2 = [cos(angle_r+pi), -sin(angle_r+pi), 0; sin(angle_r+pi), cos(angle_r+pi), 0; 0, 0, 1];
%         [F2 M2] = blade(Rr_f2*Vcg,Rr_f2*(omg+[0; 0; omg_r]),nu,pitch2_p,bladeProperties_p,rho);
%         F2 = (F2*Rr_f2);
% 
%         Rs_f1 = [cos(angle_s), -sin(angle_s), 0; sin(angle_s), cos(angle_s), 0; 0, 0, 1];
%         [F3 M3] = blade(Rs_f1*Vcg,Rs_f1*(omg+[0; 0; omg_s]),0,pitch1_d,bladeProperties_d1,rho);
%         F3 = (F3*Rs_f1);
% 
%         Rs_f2 = [cos(angle_s+pi), -sin(angle_s+pi), 0; sin(angle_s+pi), cos(angle_s+pi), 0; 0, 0, 1];
%         [F4 M4] = blade(Rs_f2*Vcg,Rs_f2*(omg+[0; 0; omg_s]),0,pitch2_d,bladeProperties_d2,rho);
%         F4 = (F4*Rs_f2);
% 
% %         T = -F1(3)-F2(3)-F3(3)-F4(3);
%         T = -F1(3)-F2(3);
%     %   Change in induced air velocity from propeller
% 
%         nd1 = -2*nu*abs(-Vcg(3) + nu)/H_p;
%         nd2 = T/(rho*area_p*H_p);
%         dNu = nd1+nd2;
% 
%     end
% 
%     %Thrust
%     Rr_f1 = [cos(angle_r), -sin(angle_r), 0; sin(angle_r), cos(angle_r), 0; 0, 0, 1];
%     [F1 M1] = blade2(Rr_f1*Vcg,Rr_f1*(omg+[0; 0; omg_r]),nu,pitch1_p,bladeProperties_p,rho);
%     F1 = (F1*Rr_f1);
%     M1 = (M1*Rr_f1);
%     
%     Rr_f2 = [cos(angle_r+pi), -sin(angle_r+pi), 0; sin(angle_r+pi), cos(angle_r+pi), 0; 0, 0, 1];
%     [F2 M2] = blade2(Rr_f2*Vcg,Rr_f2*(omg+[0; 0; omg_r]),nu,pitch2_p,bladeProperties_p,rho);
%     F2 = (F2*Rr_f2);
%     M2 = (M2*Rr_f2);
%     
%     Rs_f1 = [cos(angle_s), -sin(angle_s), 0; sin(angle_s), cos(angle_s), 0; 0, 0, 1];
%     [F3 M3] = blade2(Rs_f1*Vcg,Rs_f1*(omg+[0; 0; omg_s]),0,pitch1_d,bladeProperties_d1,rho);
%     F3 = (F3*Rs_f1);
%     M3 = (M3*Rs_f1);
%     
%     Rs_f2 = [cos(angle_s+pi), -sin(angle_s+pi), 0; sin(angle_s+pi), cos(angle_s+pi), 0; 0, 0, 1];
%     [F4 M4] = blade2(Rs_f2*Vcg,Rs_f2*(omg+[0; 0; omg_s]),0,pitch2_d,bladeProperties_d2,rho);
%     F4 = (F4*Rs_f2);
%     M4 = (M4*Rs_f2);
%     
% %     Vcg
% %     angle_r
% %     angle_s
% %     M1+M2
% %     M3+M4
% %     F1+F2
% %     F3+F4
% %     F1
% %     F2
% %     F3
% %     F4
% %     M1
% %     M2
% %     M3
% %     M4
%     
%     Fsum = F1+F2+F3+F4;
%     Msum = M1+M2+M3+M4;
% %     Mout = [Mout; F1 M1 F2 M2 F3 M3 F4 M4];
% 
%     T = -F1(3)-F2(3)-F3(3)-F4(3);
% %     T = -F1(3)-F2(3);
% %   Change in induced air velocity from propeller
% % TODO::I could use a better algorithm
%     nd1 = -2*nu*abs(-Vcg(3) + nu)/H_p;
%     nd2 = T/(rho*area_p*H_p);
%     nuDot = nd1+nd2;
% %     nu = sqrt(T/(2*rho*area_p));
% %     nuDot = 0;
%     
%     %Compute motor torque
% %     M_m = [0 0 -.03]; %(Nm)
%     M_m = [0 0 K_t_m*i_m];
%     iDot = (V - K_t_m*(omg_r-omg_s)-r_m*i_m)/L_m;
%     
%     %Force caused by gravity
%     Fg = [-m*g*sin(angles(2)), m*g*cos(angles(2))*sin(angles(1)), m*g*cos(angles(2))*cos(angles(1))];
%     
%     %CG acceleration
%     VcgDot = (1/(m))*(F1+F2+F3+F4+Fg)'-cross(omg,Vcg);
%     
%     %Rotate Inertias
%     Is_sF = Rs_f1*Is_s*Rs_f1';
%     Ir_rF = Rr_f1*Ir_r*Rr_f1';
%     
%     %Rotor rotational acceleration
%     omg_rDot = Ir_rF\((M1+M2+M_m).*[0 0 1])';
%     
%     %Stator rotational acceleration
%     omg_sDot = Is_sF\((M3+M4-M_m).*[0 0 1])';
%     
%     %CG rotational acceleration
%     omgDot = (Is_sF+Ir_rF+I_tot)\((M1+M2+M3+M4)'-cross(omg,Is_sF*(omg+[0; 0; omg_s]))-cross(omg,Ir_rF*(omg+[0; 0; omg_r]))-Is_sF*omg_sDot-Ir_rF*omg_rDot);
% 
%     %World rotational acceleration
%     anglesDot = [omg(1)+(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*tan(angles(2)); omg(2)*cos(angles(1))-omg(3)*sin(angles(1));(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*sec(angles(2))];
%     
%     %World velocity
%     Ri_f = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
%     XeDot = Ri_f*Vcg;
% 
% %     dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); omgDot(3); omg_rDot(3); omg_sDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_s; XeDot(1); XeDot(2); XeDot(3)];
%     dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); omgDot(3); omg_rDot(3); omg_sDot(3); anglesDot(1); anglesDot(2); anglesDot(3); omg_r; omg_s; XeDot(1); XeDot(2); XeDot(3); iDot];
% end