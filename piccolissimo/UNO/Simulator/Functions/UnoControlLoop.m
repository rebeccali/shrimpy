function stop = UnoControlLoop(t, X, flag) 
    global pwm_0 kpp kdp kip pwm pwm_i omg_err_prev pwm_out nu_style Vcg omg omg_r omg_b angles angle_r angle_s aoa_d_out aoa_p_out tau_p_out psi_0 kff_pitch pitch_clamp kp_pitch kd_pitch kdd_pitch target_size time F_p F_d i_m nu kpv kdv kiv omg_r_0 kdomg kpomg kiomg waypoints waypoint_num args counter x_prev v_err_prev alt_err_prev psi_err_prev alt_err_dot_prev t_prev v_i omg_i psi_i waypoint_prev z_prev report_counter v_out nu_out i_m_out thrust_p_out thrust_d_out RP_tau_p_out RP_tau_d_out M_p M_d v_clamp v v_0 K_t_m r_m pitch_p
    alpha = .5;
    if strcmpi(flag,'init')
%         v_i = 0;
        pwm_i = 0;
        omg_i = 0;
        psi_i = 0;
        t_prev = 0;
        alt_err_prev = 0;
        alt_err_dot_prev = 0;
        v_err_prev = 0;
        omg_err_prev = 0;
        psi_err_prev = 0;
        waypoint_prev = 0;
        z_prev = 0;
        x_prev = [];
        report_counter = 0;

    elseif strcmp(flag,'done')
        
    else
        %% Init normal control function
        stop = false;
        counter = counter + 1;
        report_counter = report_counter + 1;
        t_dot = t(end)-t_prev;

        %% Check if z trim is finished
        if(nnz(strcmpi(args,'trim')))
            if ~isempty(x_prev)
                if ((abs(X(4) - x_prev(4)) < .000001) && (abs(X(4)) < .0001) && (abs((X(8) - x_prev(8))/X(8))/t_dot < .0001) && (abs((X(9) - x_prev(9))/X(9))/t_dot < .0001))
                    stop = true;
                    omg_r = X(8);
                    omg_b = X(9);
                    i_m = (v - K_t_m*(omg_r-omg_b))/r_m;
%                         psidess = psidess(1:counter);
%                         v_out = v_out(1:counter);
                end
            end
        end

        %% Check if we're in a waypoint
        if(nnz(strcmpi(args,'throttle')) || nnz(strcmpi(args,'cyclic')))
            if(nnz(strcmpi(args,'spiral')))
                [pos_des, vel_des, acc_des, jrk_des] = SpiralTrajectory(t(end),0,20,[0;0;0], 0, .5, 2*pi);
            elseif(nnz(strcmpi(args,'circle')))
                pos_des = CircleTrajectory(t(end));
            elseif(sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2 + (waypoints(waypoint_num,3) - X(17))^2) < target_size)
%             disp(['At waypoint ' num2str(waypoint_num)]);
                pos_des = waypoints(waypoint_num,:);
                waypoint_num = waypoint_num + 1;
                waypoint_prev = 0;
            else
                pos_des = waypoints(waypoint_num,:);
            end
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
        v_amplitude = 0;
        psi_des = 0;
        pitch_p = [0,0];
%             if sqrt((waypoints(waypoint_num,1) - X(15))^2 + (waypoints(waypoint_num,2) - X(16))^2) > sqrt(target_size^2/2)
            if (nnz(strcmpi(args,'cyclic')) > 0 && t > 0)
                d_err = kp_pitch*(pos_des(1:2)-X(15:16));
                v_err = kd_pitch*(vel_des(1:2)-X(2:3));
                a_err = kdd_pitch*(acc_des(1:2) - [-X(11);X(10)]);
                jrk_ff = kff_pitch*jrk_des(1:2);
                tot_err = d_err + v_err + a_err + jrk_ff;
                pitch_amplitude = norm(tot_err);
                pitch_amplitude = max(-pitch_clamp,min(pitch_clamp,pitch_amplitude));
                psi_des = atan2(tot_err(2),tot_err(1)) + psi_0;

                % Change the blade angle for under actuated simulation %%%%
                pitch_p = [deg2rad(pitch_amplitude)*sin(X(12) + X(13) + X(14) + psi_des), -deg2rad(pitch_amplitude)*sin(X(12) + X(13) + X(14) + psi_des)]; % X(12) = world to flyer psi, X(13) = stator to rotor psi, X(14) flyer to stator psi
            elseif nnz(strcmpi(args,'pulsing'))
                v_amplitude = 0; %% TODO::fill me in
            end
%             end

        %% altitude control
        altitude_des = pos_des(3);
        alt_err = altitude_des - X(17,end);
        alt_err_dot = (1-alpha)*alt_err_dot_prev + (alpha)*(alt_err_prev - alt_err)/t_dot;

        omg_i = omg_i + kiomg*alt_err*t_dot;
        omg_des = omg_i + kpomg*alt_err + alt_err_dot*kdomg + omg_r_0;
        
        omg_err = omg_des-X(8);
        pwm_at_start = pwm;
        pwm_i = pwm_i + kip*(omg_err*t_dot);
        pwm = pwm_i + kpp*omg_err + (omg_err-omg_err_prev)/t_dot*kdp + pwm_0;

%         v_err = omg_des-X(8);
%         v_i = v_i + kiv*(v_err)*t_dot;
%         v = v_i + kpv*v_err + (v_err-v_err_prev)/t_dot*kdv + v_0; 

        %% mix in pulsing and clamp voltage
%         v = v + v_amplitude*sign(cos(X(14)+X(12)+psi_des));
%         if v > v_clamp
%             v = v_clamp;
%         elseif v < -v_clamp
%             v = -v_clamp;
%         end

        if pwm > 1
            pwm = 1;
        elseif pwm < -1
            pwm = -1;
        end

        % Store current values for next calculation
        t_prev = t(end);
        x_prev = X; 
        alt_err_dot_prev = alt_err_dot;
        alt_err_prev = alt_err;
%         v_err_prev = v_err;
        omg_err_prev = omg_err;
        waypoint_prev = 1;

        % Store any desired variables
        v_out(counter) = v;
        nu_out(:,counter) = nu;
        i_m_out(counter) = i_m;
        thrust_p_out(counter) = F_p(3);
        thrust_d_out(counter) = F_d(3);
        RP_tau_p_out(counter) = sqrt(M_p(1)*M_p(1) + M_p(2)*M_p(2));
        RP_tau_d_out(counter) = sqrt(M_d(1)*M_d(1) + M_d(2)*M_d(2));
        tau_p_out(:,counter) = M_p;
        pwm_out(counter) = pwm_at_start;
        
        %% Stuff to help find AOA
            Vcg = X(2:4);
            omg = X(5:7);
            omg_r = X(8);
            omg_b = X(9);
            angles = X(10:12);
            angle_r = X(13);
            angle_s = X(14);
            if(nu_style == 1)
                nu = NumericalFixedPointInflow(nu);
            elseif(nu_style == 0)
                nu = X(1);
            end
            [~, ~, ~, ~, aoa_p, aoa_d] = ComputeAero();
            aoa_p_out(:,:,end+1) = aoa_p;
            aoa_d_out(:,:,end+1) = aoa_d;

        %%Display Simulation Status
        if toc > 1
            figure(9);
            if((exist('pitch_amplitude', 'var') && exist('pitch_clamp', 'var'))&&(pitch_amplitude == pitch_clamp || pitch_amplitude == -pitch_clamp))
                plot_style = '.r';
            else
                plot_style = '.b';
            end
            plot3(X(15),X(16),X(17),plot_style,'Clipping','off','linewidth',3);
            plot3(pos_des(1),pos_des(2),pos_des(3),'.g','Clipping','off','linewidth',3);
            drawnow();

            disp(['Sim time: ',num2str(t(end),5),', time per calc: ',num2str(toc/report_counter,5), ', time left: ' num2str((time-t(end))*toc/report_counter*1000)]);
            tic
            report_counter = 0;
        end
    end
end