function [t,X] = simulatePiccolissimoV2(varargin)
%simulatePiccolissimoV2 simulates a piccolissimo with a vertical, offset
%motor.  It can simulate with/without a synchronous driver and with/without
%pulsing.
% 
% X = [phi_b, theta_b, psi_b, psi_p, phi_b_dot, theta_b_dot, psi_b_dot, psi_p_dot];

pulsing = getArgVal('pulsing', varargin, true);
synchronous = getArgVal('synchronous', varargin, true);
mass = getArgVal('mass', varargin, .0042);
t_max = getArgVal('t_max',varargin, 10);
plot_hold = getArgVal('plot_hold', varargin, false);

k_t_prop = 2.66e-11*(60/(2*pi))^2;%2.33e-11*(60/(2*pi))^2; %N/RPM^2 -> N/(rad/s)^2
k_d_prop = 1.15e-13*(60/(2*pi))^2;%1.98e-13*(60/(2*pi))^2; %Nm/RPM^2 -> Nm/(rad/s)^2
k_d_body = 2.15e-09;%3.0245e-09; %Nm/(rad/s)^2
inertia_prop = 2.7255e-9; %kg*m*m
inertia_rotor = .00025*.0025^2; %.25g, 5mm dia
I_p = [.00025/12*(3*.0025^2+.009),0,0; 0,.00025/12*(3*.0025^2+.009),0; 0,0,inertia_prop+inertia_rotor];
I_b = [410,0,0;0,550,0;0,0,820]*1e-9;
v_batt = 4.2;
k_v = 20500; % 20000 = cheerson motor
R = 1.97; % 1.97 = cheerson motor
motor_offset = .0046; % motor location along y direction 
psi_offset = -pi/4;

k_e = 1/(k_v/60*2*pi);
hover_thrust = mass*9.81;
omega_p0 = sqrt(hover_thrust/k_t_prop);
tau_m_base = k_d_prop*omega_p0*omega_p0;
omega_b0 = sqrt(tau_m_base/k_d_body); %rad/s
v_base = (omega_p0+omega_b0)*k_e + tau_m_base/k_e*R;

if(synchronous)
    pwm_base = v_base/v_batt;
else
    pwm_base = (v_base - (omega_p0+omega_b0)*k_e)/(v_batt-(omega_p0+omega_b0)*k_e);
end

x0 = [0,0,0,0,0,0,omega_b0,-omega_p0];

t_out = [];
tau_offset_out = [];
tau_gyro_out = [];

options = odeset('OutputFcn',@simOutputFun,'RelTol',1e-5,'AbsTol',1e-8);
[t,X] = ode113(@simFun,0:0.001:t_max,x0,options);

plotSim(t,X,plot_hold);
x_ss_rate = mean(tau_offset_out(1,:))/(I_b(3,3)*mean(X(:,7)) + I_p(3,3)*mean(X(:,8)));
y_ss_rate = -mean(tau_offset_out(2,:))/(I_b(3,3)*mean(X(:,7))+ I_p(3,3)*mean(X(:,8)));
xy_ss_rate = sqrt(mean(tau_offset_out(1,:)).^2 + mean(tau_offset_out(2,:)).^2)/(I_b(3,3)*mean(X(:,7))+ I_p(3,3)*mean(X(:,8)));
figure(5); plot(t_out,tau_offset_out,t_out,tau_gyro_out);
figure(1); hold all; 
plot([0,t_out(end)],[0,t_out(end)*x_ss_rate]);
plot([0,t_out(end)],[0,t_out(end)*y_ss_rate]);
plot([0,t_out(end)],[0,t_out(end)*xy_ss_rate]);
figure(4); hold all;
plot([0,t_out(end)],[x_ss_rate, x_ss_rate],[0,t_out(end)],[y_ss_rate, y_ss_rate],[0,t_out(end)],[xy_ss_rate, xy_ss_rate])

mean(tau_offset_out(1,:))
mean(tau_gyro_out(1,:))
mean(tau_offset_out(1,:)-tau_gyro_out(1,:))

mean(tau_offset_out(2,:))
mean(tau_gyro_out(2,:))
mean(tau_offset_out(2,:)-tau_gyro_out(2,:))


    function x_dot = simFun(t,x)
        phi_b = x(1);
        theta_b = x(2);
        psi_b = x(3);
        psi_p = x(4);
        phi_b_dot = x(5);
        theta_b_dot = x(6);
        psi_b_dot = x(7);
        psi_p_dot = x(8);

        %% Compute motor pwm
        if(pulsing)
            pwm_now = wrapTo2Pi(psi_b + psi_offset) > 2*pi*(1-pwm_base);
        else %not pulsing
            pwm_now = pwm_base;
        end

        %% Compute motor voltage
        if(synchronous)
            v_now = v_batt*pwm_now;
        else % asynchrnous
            v_now = v_batt*pwm_now + (1-pwm_now)*(psi_b_dot-psi_p_dot)*k_e;
        end
        i_now = (v_now - k_e*(psi_b_dot-psi_p_dot))/R; % inductance is too fast

        %% Compute torques
        [tau_m, tau_b_a, tau_p_a, tau_b_gyro, tau_p_gyro, tau_offset] = computeTorques(i_now, x);
        tau_offset_vect = [tau_offset;0;0];
        [Rb_f, Rp_f] = computeRotations(x);

        %% Compute accelerations
        psi_p_dot_dot = (tau_p_a-tau_m)/(I_p(3,3));
        psi_b_dot_dot = (tau_b_a+tau_m)/(I_b(3,3));
        ptp_b_dot_dot = (Rp_f*I_p*Rp_f'+Rb_f*I_b*Rb_f')\(Rb_f'*tau_offset_vect-tau_b_gyro-tau_p_gyro);

        %% Populate state vector
        x_dot = [phi_b_dot; theta_b_dot; psi_b_dot; psi_p_dot; ...
            ptp_b_dot_dot(1); ptp_b_dot_dot(2); psi_b_dot_dot; psi_p_dot_dot];
    end

    function status = simOutputFun(t,x,flag)
        if isempty(flag)
            if(abs(x(1) > pi/2 || abs(x(2)) > pi/2))
                status = 1;
                return;
            end
            
            %% Compute torques
            [tau_m, tau_b_a, tau_p_a, tau_b_gyro, tau_p_gyro, tau_offset] = computeTorques(nan, x);
            tau_offset_vect = [tau_offset;0;0];
            [Rb_f, Rp_f] = computeRotations(x);
            
            %% Store stuff
            t_out(end+1) = t;
            tau_offset_out(:,end+1) = Rb_f'*tau_offset_vect;
            tau_gyro_out(:,end+1) = tau_b_gyro + tau_p_gyro;
        end
        status = 0;
    end

    function plotSim(t,X,plot_hold)
        figure(1);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,1),t,X(:,2),t,sqrt(X(:,1).*X(:,1)+X(:,2).*X(:,2)));
        xlabel('Time (s)');
        ylabel('Angle (rad)');
        title('Body Angles vs. Time');
        
        figure(2);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,7));
        xlabel('Time (s)');
        ylabel('Body Z rotation rate (rad/s)');
        title('Body Rate vs. Time');
        
        figure(3);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,8));
        xlabel('Time (s)');
        ylabel('Propeller Z rotation rate (rad/s)');
        title('Propeller Rate vs. Time');
        
        figure(4);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,5),t,X(:,6),t,sqrt(X(:,5).^2 + X(:,6).^2));
        xlabel('Time (s)');
        ylabel('Angular Velocity (rad/s)');
        title('Body Rates vs. Time');
    end

    function [Rb_f, Rp_f] = computeRotations(x)
        psi_b = x(3);
        psi_p = x(4);
        
        Rb_f = [cos(psi_b), -sin(psi_b), 0; sin(psi_b), cos(psi_b), 0; 0, 0, 1];
        Rp_f = [cos(psi_p), -sin(psi_p), 0; sin(psi_p), cos(psi_p), 0; 0, 0, 1];
%         Rb_f = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%         Rp_f = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    end

    function [tau_m, tau_b_a, tau_p_a, tau_b_gyro, tau_p_gyro, tau_offset] = computeTorques(i_now, x)
        phi_b_dot = x(5);
        theta_b_dot = x(6);
        psi_b_dot = x(7);
        psi_p_dot = x(8);
        
        tau_m = i_now*k_e;
        tau_b_a = -k_d_body*psi_b_dot*psi_b_dot;
        tau_p_a = k_d_prop*psi_p_dot*psi_p_dot;
        tau_b_gyro = cross3([phi_b_dot;theta_b_dot;0],I_b*[phi_b_dot;theta_b_dot;psi_b_dot]);
        tau_p_gyro = cross3([phi_b_dot;theta_b_dot;0],I_p*[phi_b_dot;theta_b_dot;psi_p_dot]);
        tau_offset = k_t_prop*psi_p_dot*psi_p_dot*motor_offset;
    end
end