function [t,X] = simulatePiccolissimoV1(state_transition, varargin)
%simulatePiccolissimoV2 simulates a piccolissimo with an angled
%motor.  It can simulate with/without a synchronous driver and with/without
%pulsing.
% 
% X = [u, v, p, q, phi, theta, r_p, r_b, psi_p, psi_b, x, y]
% Piccolissimo_V7:
% state_transition = [ -1.1204 -3.0867e-08 -0.0044094 0.0041307 1.0897e-14 -9.8098; 5.2105e-06 -1.1204 -0.0041305 -0.0044088 9.8098 -5.1417e-15; -86.3872 -40.5605 -1.7128 455.9048 7.6265e-13 7.516e-13; 40.5608 -86.3871 -455.9048 -1.7126 6.0061e-13 6.2651e-13; 0 0 1 0 0 0; 0 0 0 1 0 0]
% state_transition = [ -1.1204 -3.0867e-08 -0.0044094 0.0041307 1.0897e-14 -9.8098; 5.2105e-06 -1.1204 -0.0041305 -0.0044088 9.8098 -5.1417e-15; -86.3872 0 -1.7128 455.9048 7.6265e-13 7.516e-13; 0 -86.3871 -455.9048 -1.7126 6.0061e-13 6.2651e-13; 0 0 1 0 0 0; 0 0 0 1 0 0]
% setup_flyer_piccolissimo_V7 %-4.5 deg beta_p %.5 deg beta_b
% state_transition = [-1.10009647804938,-2.58301125878470e-05,-0.00341156108488691,0.00401831247791120,1.05443811380609e-13,-9.80999836500000;2.48836865896714e-05,-1.10009883008831,-0.00400075649448386,-0.00340334227113705,9.80999836500004,-5.55684448873283e-14;-82.1701695805877,-39.8275059934519,-1.56814517143908,456.795144377756,6.89998204475747e-12,6.78081458150171e-12;39.8273883282910,-82.1701428837520,-456.797906286225,-1.56898598571314,5.07070489079015e-12,4.95851670515334e-12;0,0,1.00000000000000,0,0,0;0,0,0,1.00000000000000,0,0];
% setup_flyer_piccolissimo_V7 %-4.5 deg beta_p %.5 deg beta_b % No c
% state_transition = [-1.10009647804938,-2.58301125878470e-05,-0.00341156108488691,0.00401831247791120,1.05443811380609e-13,-9.80999836500000;2.48836865896714e-05,-1.10009883008831,-0.00400075649448386,-0.00340334227113705,9.80999836500004,-5.55684448873283e-14;-82.1701695805877,0,-1.56814517143908,456.795144377756,6.89998204475747e-12,6.78081458150171e-12;0,-82.1701428837520,-456.797906286225,-1.56898598571314,5.07070489079015e-12,4.95851670515334e-12;0,0,1.00000000000000,0,0,0;0,0,0,1.00000000000000,0,0];
%
% Piccolissimo_V11
% state_transition = [-1.40212373746823,-1.75202605714352e-06,-0.00570824087720790,-0.000496773009854184,5.53273807526393e-14,-9.80999836500004;9.29033512644659e-07,-1.40211228446611,0.000500633508900422,-0.00570589607940194,9.80999836500009,-2.98857734633648e-14;-86.3338489040006,2.19799865165779,-1.08886006588965,416.023382161326,1.90344021906833e-12,2.04430080077836e-12;-2.19813550951060,-86.3332950241713,-416.024855091881,-1.08914264206303,1.23786269734535e-11,1.23495829994388e-11;0,0,1.00000000000000,0,0,0;0,0,0,1.00000000000000,0,0];
% Piccolissimo_V11 paper
% state_transition = [-1.33783549862383,3.60082297200755e-05,-0.00448341311616154,-0.000264760781081369,4.43673144444092e-14,-9.80999836500006;-2.92479923190857e-05,-1.33784241483851,0.000292327179283993,-0.00446310298707323,9.80999836500009,-2.84185147347971e-14;-21.9535440231964,1.68995378921132,-1.10710841013012,453.820282775770,1.81690909498165e-13,2.24371166038912e-13;-1.69091802966630,-21.9523111640813,-453.824853996392,-1.11159149286849,4.54704893901275e-12,4.61661470116886e-12;0,0,1.00000000000000,0,0,0;0,0,0,1.00000000000000,0,0];
if nargin == 0
    state_transition = [-1.33783549862383,3.60082297200755e-05,-0.00448341311616154,-0.000264760781081369,4.43673144444092e-14,-9.80999836500006;-2.92479923190857e-05,-1.33784241483851,0.000292327179283993,-0.00446310298707323,9.80999836500009,-2.84185147347971e-14;-21.9535440231964,1.68995378921132,-1.10710841013012,453.820282775770,1.81690909498165e-13,2.24371166038912e-13;-1.69091802966630,-21.9523111640813,-453.824853996392,-1.11159149286849,4.54704893901275e-12,4.61661470116886e-12;0,0,1.00000000000000,0,0,0;0,0,0,1.00000000000000,0,0];
end

pulsing = getArgVal('pulsing', varargin, true);
synchronous = getArgVal('synchronous', varargin, true);
mass = getArgVal('mass', varargin, .00447);
t_max = getArgVal('t_max',varargin, 10);
plot_hold = getArgVal('plot_hold', varargin, false);
motor_angle = getArgVal('motor_angle', varargin, deg2rad(0));
motor_offset = getArgVal('motor_offset', varargin, .00807);
k_t_prop_mult = getArgVal('k_t_prop_mult', varargin, 1);
x0 = getArgVal('x0',varargin, [0,0,0,0,0,0,0,0,0,0,0,0]);

k_t_prop = 2.66e-11*(60/(2*pi))^2*k_t_prop_mult;%2.33e-11*(60/(2*pi))^2; %N/RPM^2 -> N/(rad/s)^2
k_d_prop = 1.15e-13*(60/(2*pi))^2*k_t_prop_mult;%1.98e-13*(60/(2*pi))^2; %Nm/RPM^2 -> Nm/(rad/s)^2
% k_t_prop = 2*2.66e-11*(60/(2*pi))^2;%FAKE
% k_d_prop = 2*1.15e-13*(60/(2*pi))^2;%FAKE
k_d_body = 2.15e-09;%3.0245e-09; %Nm/(rad/s)^2
% inertia_prop = 2.7255e-9; %kg*m*m
% inertia_rotor = .00025*.0025^2; %.25g, 5mm dia
% I_p = [.00025/12*(3*.0025^2+.009),0,0; 0,.00025/12*(3*.0025^2+.009),0; 0,0,inertia_prop+inertia_rotor];
% I_b = [410,0,0;0,550,0;0,0,820]*1e-9;
I_p = [ .1+1.36,0,0;0,2.71+1.36,0;0,0, 2.73+0.97]*1e-9; %Cheerson prop plus motor rotor %Rotor inertia at rotor cg in flyer frame 
% I_b = [463,0,0;0,468,0;0,0,783]*1e-9; %Stator inertia at stator cg in flyer frame
I_b = 2.5*[383, 0, 0; 0, 439, 0; 0, 0, 697]*1e-9;
v_batt = 4.2;
k_v = 20500; % 20000 = cheerson motor
R = 1.97; % 1.97 = cheerson motor
psi_offset = pi/4+pi+pi/16;%pi+pi/4-pi/32+pi/128;%-pi/2+1/8*pi+1/128*pi;
psi_angle = 0;
Ra_b = [cos(psi_angle), -sin(psi_angle),0; sin(psi_angle), cos(psi_angle),0;0,0,1];
R_batt = 2.0; %(4.2-2.7)/1.25; % dropping to 2.7V, maybe pulling 1.25A when pulsing!!!

k_e = 1/(k_v/60*2*pi);
hover_thrust = mass*9.81;
omega_p0 = sqrt(hover_thrust/k_t_prop);
tau_m_base = k_d_prop*omega_p0*omega_p0;
omega_b0 = sqrt(tau_m_base/k_d_body); %rad/s
i_now = tau_m_base/k_e;
v_base = (omega_p0+omega_b0)*k_e + i_now*R; % average voltage across motor
v_now = v_base;
pwm_sync2 = roots([i_now*R_batt, -v_batt, v_base]);
pwm_sync = pwm_sync2(pwm_sync2 >=0 & pwm_sync2 <=1);
if(isempty(pwm_sync))
    error('PWM out of limits');
end
v_now_batt = v_batt-i_now*pwm_sync*R_batt;

if(synchronous)
    pwm_base = v_base/v_now_batt; % should equal pwm_sync
else
    pwm_base = (v_base - (omega_p0+omega_b0)*k_e)/(v_now_batt-(omega_p0+omega_b0)*k_e);
end

% x0 = [0,0,0,0,0,0,omega_p0,-omega_b0,0,0,0,0];
x0(7) = omega_p0;
x0(8) = -omega_b0;

t_out = [];
tau_offset_out = [];
v_out = [];
i_out = [];
v_batt_out = [];

options = odeset('OutputFcn',@simOutputFun,'MaxStep',.0001,'RelTol',1e-5,'AbsTol',1e-8);
% options = odeset('AbsTol',1e-3,'RelTol',1e-2,'MaxStep',.0001,'OutputFcn',@simOutputFun);
[t,X] = ode113(@simFun,0:0.0001:t_max,x0,options);

plotSim(t,X,plot_hold);

%% Additional plot stuff
tau_x = mean(tau_offset_out(1,:));
tau_y = mean(tau_offset_out(2,:));
tau_xy = sqrt(tau_x.^2 + tau_y.^2);
I_b_xy = mean([I_b(1,1) I_b(2,2)]);
I_p_xy = mean([I_p(1,1) I_p(2,2)]);
r_b = mean(X(:,8));
r_p = mean(X(:,7));
p_ss = -tau_y/(I_b(3,3)*r_b + I_p(3,3)*r_p);
q_ss = tau_x/(I_b(3,3)*r_b+ I_p(3,3)*r_p);
pq_ss = -tau_xy/(I_b(3,3)*r_b+ I_p(3,3)*r_p);
k_u = -1/(I_b_xy + I_p_xy)*tau_x/state_transition(3,1);
k_v = -1/(I_b_xy + I_p_xy)*tau_y/state_transition(3,1);
k_uv = -1/(I_b_xy + I_p_xy)*tau_xy/state_transition(3,1);
k_phi = -state_transition(1,1)*k_v/state_transition(2,5);
k_theta = state_transition(1,1)*k_u/state_transition(2,5);
k_phitheta = -state_transition(1,1)*k_uv/state_transition(2,5);
t_x_ss = k_phi/p_ss;
t_y_ss = k_theta/q_ss;
t_xy_ss = k_phitheta/pq_ss;

figure(2); 
hold all;
plot([0,t_x_ss,t_out(end)],[0,k_phi,k_phi]);
plot([0,t_y_ss,t_out(end)],[0,k_theta,k_theta]);
plot([0,t_xy_ss,t_out(end)],[0,k_phitheta,k_phitheta]);

figure(33); 
hold all;
plot([0,t_out(end)],[k_u,k_u]);
plot([0,t_out(end)],[k_v,k_v]);
plot([0,t_out(end)],[k_uv,k_uv]);

figure(37);
if(plot_hold)
    hold all;
else
    clf;
end
plot(1:size(tau_offset_out,2),tau_offset_out)
hold all;
% plot(t_out,tau_offset_out);
plot([0,size(tau_offset_out,2)],[mean(tau_offset_out,2), mean(tau_offset_out,2)]);
plot(1:size(tau_offset_out,2),sqrt(tau_offset_out(1,:).*tau_offset_out(1,:)+tau_offset_out(2,:).*tau_offset_out(2,:)));
xlabel('Time (timestep)');
ylabel('Torque (Nm)');
title('Torque vs. Timestep');

figure(38);
if(plot_hold)
    hold all;
else
    clf;
end
plot(t_out,v_out,t_out,v_batt_out,t_out,i_out);
legend('V_{out}', 'V_{batt}', 'I');
xlabel('Time (s)');
ylabel('Voltage (Volts)');
title('Voltage vs. Time');

    function x_dot = simFun(t,x)
        u = x(1);
        v = x(2);
        p = x(3);
        q = x(4);
        phi = x(5);
        theta = x(6);
        r_p = x(7);
        r_b = x(8);
        psi_p = x(9);
        psi_b = x(10);

        %% Compute motor pwm
        if(pulsing)
            pwm_now = wrapTo2Pi(psi_b + psi_offset) < 2*pi*pwm_base;
        else %not pulsing
            pwm_now = pwm_base;
        end
        i_now = (v_batt*pwm_now-k_e*(-r_b+r_p))/(pwm_now*pwm_now*R_batt+R);
        v_now_batt = v_batt-R_batt*i_now*pwm_now;

        %% Compute motor voltage
        if(synchronous)
%             v_now = v_batt*pwm_now; % without battery R loss
            v_now = pwm_now*v_now_batt;
        else % asynchrnous
            v_now = pwm_now*v_now_batt + (1-pwm_now)*(-r_b+r_p)*k_e;
        end
%         i_now = (v_now - k_e*(-r_b+r_p))/(R + R_batt*pwm_now); % inductance is too fast

        %% Compute torques
        [tau_m, tau_b_a, tau_p_a, tau_offset] = computeTorques(i_now, x);
        tau_offset_vect = [tau_offset;0;0];
        [Rb_f, Rp_f] = computeRotations(x);
        offset_sensitivity = (Rp_f*I_p*Rp_f'+Rb_f*I_b*Rb_f')\(Rb_f*tau_offset_vect);
        
        %% Compute forces
        F_angle = Rb_f'*Ra_b'*[0;mass*9.81*sin(motor_angle);0];

        %% Compute accelerations
        r_p_dot = (tau_p_a+tau_m)/(I_p(3,3));
        r_b_dot = (tau_b_a-tau_m)/(I_b(3,3));
        Bu = [F_angle(1)/mass; F_angle(2)/mass; offset_sensitivity(1); offset_sensitivity(2); 0; 0];

        %% Populate state vector
        x_dot = [(state_transition*x(1:6) + Bu); r_p_dot; r_b_dot; r_p; r_b; u; v];
    end

    function status = simOutputFun(t,x,flag)
        if isempty(flag)
            if(abs(x(1) > pi/2 || abs(x(2)) > pi/2))
                status = 1;
                return;
            end
            
            %% Compute torques
            [tau_m, tau_b_a, tau_p_a, tau_offset] = computeTorques(nan, x);
            tau_offset_vect = [tau_offset;0;0];
            [Rb_f, Rp_f] = computeRotations(x);
            
            %% Store stuff
            t_out = [t_out, t(end)];
            tau_offset_out = [tau_offset_out, Rb_f*tau_offset_vect];
            
            v_out = [v_out, v_now];
            i_out = [i_out, i_now];
            v_batt_out = [v_batt_out, v_now_batt];
        end
        status = 0;
    end

    function plotSim(t,X,plot_hold)
        figure(2);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,5),t,X(:,6),t,sqrt(X(:,5).*X(:,5)+X(:,6).*X(:,6)));
        xlabel('Time (s)');
        ylabel('Angle (rad)');
        title('Body Angles vs. Time');
        
        figure(32);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,8));
        xlabel('Time (s)');
        ylabel('Body Z rotation rate (rad/s)');
        title('Body Rate vs. Time');
        
%         figure(33);
%         if(plot_hold)
%             hold all;
%         else
%             clf;
%         end
%         plot(t,X(:,7));
%         xlabel('Time (s)');
%         ylabel('Propeller Z rotation rate (rad/s)');
%         title('Propeller Rate vs. Time');
        
        figure(34);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,3),t,X(:,4),t,sqrt(X(:,3).^2 + X(:,4).^2));
        xlabel('Time (s)');
        ylabel('Angular Velocity (rad/s)');
        title('Body Rates vs. Time');
        
        figure(33);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,1),t,X(:,2),t,sqrt(X(:,1).^2 + X(:,2).^2));
        xlabel('Time (s)');
        ylabel('Velocity (m/s)');
        title('Velocity vs. Time');
        
        figure(1);
        if(plot_hold)
            hold all;
        else
            clf;
        end
        plot(t,X(:,11),t,X(:,12));
        xlabel('Time (s)');
        ylabel('Position (m)');
        title('Position vs. Time');
    end

    function [Rb_f, Rp_f] = computeRotations(x)
        psi_p = x(9);
        psi_b = x(10);
        
        Rb_f = [cos(psi_b), -sin(psi_b), 0; sin(psi_b), cos(psi_b), 0; 0, 0, 1];
        Rp_f = [cos(psi_p), -sin(psi_p), 0; sin(psi_p), cos(psi_p), 0; 0, 0, 1];
%         Rb_f = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%         Rp_f = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    end

    function [tau_m, tau_b_a, tau_p_a, tau_offset] = computeTorques(i_now, x)
        r_p = x(7);
        r_b = x(8);
        
        tau_m = i_now*k_e;
        tau_b_a = k_d_body*r_b*r_b;
        tau_p_a = -k_d_prop*r_p*r_p;
        tau_offset = -k_t_prop*r_p*r_p*motor_offset;
    end
end