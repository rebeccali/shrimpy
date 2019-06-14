function stop = FlyerControlLoop(t, X, flag) 
    global target_size time F_p F_d i_m nu kpv kdv kiv omg0 kdomg kpomg kiomg waypoints waypoint_num args counter xOld errvPrev erraltPrev errpsiPrev erraltDotPrev tPrev V_i omg_i psi_i waypoint_old ZPrev reportCounter Vs Nu_out i_m_out T_p_out T_d_out RP_tau_p_out RP_tau_d_out M_p M_d V_clamp V V0 K_t_m r_m pitch_p
    alpha = .5;
    if strcmpi(flag,'init')
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
                if ((abs(X(4) - xOld(4)) < .000001) && (abs(X(4)) < .0001) && (abs((X(8) - xOld(8))/X(8))/dt < .0001) && (abs((X(9) - xOld(9))/X(9))/dt < .0001))
                    stop = true;
                    omg_r = X(8);
                    omg_b = X(9);
                    i_m = (V - K_t_m*(omg_r-omg_b))/r_m;
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
%                   amplitude = min(kpcyc*sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2),30/256*3.7); % proportional pulsing amplitude up to 1
%                     amplitude = .37/2 * min(1,(kpcyc*sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2)) > target_size);% 30/256*4.2
%                 amplitude = .37/2;
%                     amplitude = V_clamp + V0; %because V0 is negative
                angleToTarget = atan2(waypoints(waypoint_num,2) - X(16),waypoints(waypoint_num,1) - X(15)); % calculate angle to desired waypoint in world
%                     psides = angleToTarget + psi0;

                % for piccolissimo pulsing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 angleCurrent = -atan2(X(3),X(2)) - mod(X(12),2*pi); % velocity angle in flyer - world to flyer angle = velocity angle in world
%                 angleCurrent = mod(angleCurrent,2*pi); 
%                 if angleCurrent > pi
%                     angleCurrent = angleCurrent - 2*pi;
%                 end
%                 angleError = mod(angleToTarget-angleCurrent,2*pi);
%                 if angleError > pi
%                     angleError = angleError - 2*pi;
%                 end
%                 % max speed .2 m/s
%                 psides = angleToTarget+0.1*angleError + psi0;
%                 psides = psi0;
% 
%                 psidess(counter) = psides;

                % Change the blade angle for under actuated simulation %%%%
                    pitch_p = [deg2rad(amplitude)*sin(X(12) + X(13) + psides), -deg2rad(amplitude)*sin(X(12) + X(13) + psides)]; % X(12) = pilot frame psi, X(13) = rotor psi
            end

%             end

        %% altitude control
        altitude_des = waypoints(waypoint_num,3);
        erralt = altitude_des - X(17,end);
        erraltdot = (1-alpha)*erraltDotPrev + (alpha)*(erraltPrev - erralt)/dt;

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
        Nu_out(:,counter) = nu;
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