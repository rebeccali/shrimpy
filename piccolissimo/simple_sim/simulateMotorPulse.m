function [ t,X, t_out, i_out, v_out, tau_out, thrust_out] = simulateMotorPulse(varargin)
%simulateMotorPulse simulates the angled pulsing motor of piccolissimo
%   x = [dist_x, dist_y, velocity_x, velocity_y, omega]

pulsing = getArgVal('pulsing', varargin, true);
synchronous = getArgVal('synchronous', varargin, true);
mass = getArgVal('mass', varargin, .0042);
t_max = getArgVal('t_max',varargin, 1);
plot_hold = getArgVal('plot_hold', varargin, false);

pulse_start = 0;
psi_offset = -pi/4;

k_t_prop = 2.66e-11*(60/(2*pi))^2;%2.33e-11*(60/(2*pi))^2; %N/RPM^2 -> N/(rad/s)^2
k_d_prop = 1.15e-13*(60/(2*pi))^2;%1.98e-13*(60/(2*pi))^2; %Nm/RPM^2 -> Nm/(rad/s)^2
k_d_body = 2.15e-09;%3.0245e-09; %Nm/(rad/s)^2
inertia_prop = 2.7255e-9; %kg*m*m
inertia_rotor = .00025*.0025^2; %.25g, 5mm dia
inertia = inertia_prop+inertia_rotor;
motor_angle = deg2rad(8);
vehicle_mass = mass; %kg
%tau_m_amp = 30/200;
v_batt = 4.2;
pwm_amp = 30/256;
coeff_drag_body = 1; % flat plate = 1.28, sphere = .07 to .5
area_body = .039*.018; % small .0342*.0177, large .039*.018
rho = 1.2041;
k_v = 20500; % 20000 = cheerson motor
R = 1.97; % 1.97 = cheerson motor

body_theta = 0;
k_e = 1/(k_v/60*2*pi);
hover_thrust = vehicle_mass*9.81/cos(motor_angle);
omega_0 = sqrt(hover_thrust/k_t_prop);
tau_m_base = k_d_prop*omega_0*omega_0;
tau_m = tau_m_base;
body_rate = sqrt(tau_m_base/k_d_body) %rad/s
v_base = (omega_0+body_rate)*k_e + tau_m_base/k_e*R;

if(synchronous)
    pwm_base = v_base/v_batt;
else
    pwm_base = (v_base - (omega_0+body_rate)*k_e)/(v_batt-(omega_0+body_rate)*k_e);
end

x0 = [0, 0, 0, 0, omega_0];

t_out = [];
i_out = [];
v_out = [];
tau_out = [];
thrust_out = [];

options = odeset('OutputFcn',@simOutputFun,'RelTol',1e-5,'AbsTol',1e-8);
[t,X] = ode113(@simFun,0:.001:t_max,x0,options);

    function x_dot = simFun(t,x)
        body_theta = t(end)*body_rate;
        if(t(end) > pulse_start)
%             v_applied = min(pwm_base*v_batt + sign(sin(body_theta))*pwm_amp*v_batt,v_batt);

%             pwm_current = min(pwm_base + sign(sin(body_theta))*pwm_amp,1);
%             v_applied = pwm_current*v_batt + (1-pwm_current)*(x(5)+body_rate)*k_e;
            
            %tau_m = tau_m_base + sign(sin(body_theta))*tau_m_base*tau_m_amp;
            
            if(pulsing)
                pwm_now = wrapTo2Pi(body_theta + psi_offset) > 2*pi*(1-pwm_base);
            else %not pulsing
                pwm_now = pwm_base;
            end
        else
%             v_applied = min(v_base,v_batt);
            %tau_m = tau_m_base;
            pwm_now = pwm_base;
        end
        
        %% Compute motor voltage
        if(synchronous)
            v_now = v_batt*pwm_now;
        else % asynchrnous
            v_now = v_batt*pwm_now + (1-pwm_now)*(x(5)+body_rate)*k_e;
        end
        
        v_out = [v_out; v_now];
        i_now = (v_now - k_e*(x(5)+body_rate))/R; % inductance is too fast
        t_out = [t_out; t];
        i_out = [i_out; i_now];
        thrust_out = [thrust_out; x(5)*x(5)*k_t_prop];
        tau_m = i_now*k_e;
        tau_out = [tau_out, tau_m]; %  - k_d_prop*x(5)*x(5)
        [acceleration_x, acceleration_y] = forceOde(x(3),x(4),x(5));
        x_dot = [x(3); x(4); acceleration_x; acceleration_y; torqueOde(tau_m,x(5))];
    end

    function status = simOutputFun(t,x,flag)
        if isempty(flag)
%             if sqrt(x(3)*x(3)+x(4)*x(4)) > .04*3 % if traveling at 3 body lengths per second
%                 status = 1;
%                 return;
%             end
                
        end
        status = 0;
    end
        
    function omega_dot = torqueOde(tau, omega_in)
        tau_eff = tau - k_d_prop*omega_in*omega_in;
        omega_dot = tau_eff/inertia;
    end

    function [acceleration_x, acceleration_y] = forceOde(velocity_x,velocity_y,omega)
        acceleration_x = (sin(motor_angle)*sin(body_theta) * k_t_prop*omega*omega - sign(velocity_x)*.5*rho*velocity_x*velocity_x*coeff_drag_body*area_body)/vehicle_mass;
        acceleration_y = (sin(motor_angle)*cos(body_theta) * k_t_prop*omega*omega - sign(velocity_y)*.5*rho*velocity_y*velocity_y*coeff_drag_body*area_body)/vehicle_mass;
    end

figure(1);
plot(t,X(:,1),t,X(:,2));
title('Positions');
figure(2);
plot(t,X(:,3),t,X(:,4));
title('Velocities');
figure(3);
plot(t,X(:,5));
title('Omega');
figure(4);
plot(t_out,i_out,t_out,v_out,[t_out(1) t_out(end)],[v_batt v_batt]);
title('Current and Voltage');
figure(5);
plot(t_out,tau_out);
title('Torque');
figure(6);
plot(t,X(:,5).*X(:,5)*k_t_prop);
title('Thrust');

%% Plost FFTs
    % Force
    figure(7);
%     hold all;
    [freq, amp] = trialFFT(t,X(:,5).*X(:,5)*k_t_prop);
    plot(freq,amp);
    title('Force FFT');
    xlim([0 100]);
    
    % Force X
    figure(8);
%     hold all;
    [freq, amp] = trialFFT(t_out,tau_out'*1000);
    plot(freq,amp);
    title('Torque FFT');
    xlim([0 100]);
    
    % Force X
    figure(9);
%     hold all;
    [freq, amp] = trialFFT(t,X(:,5)/2/pi*60);
    plot(freq,amp);
    title('Speeds FFT');
    xlim([0 100]);
    
    % Voltage
    figure(10);
%     hold all;
    [freq, amp] = trialFFT(t_out,v_out);
    plot(freq,amp);
    title('Voltage FFT');
    xlim([0 100]);
    
    % Current
    figure(11);
%     hold all;
    [freq, amp] = trialFFT(t_out,i_out);
    plot(freq,amp);
    title('Current FFT');
    xlim([0 100]);


end

function [freq, amp] = trialFFT(time,data)
    [time, sort_inds] = sort(time);
    data = data(sort_inds);
    [time,unique_inds,~] = unique(time);
    data = data(unique_inds);
    time_interp = 0:.00005:time(end);
    data_interp = interp1q(time,data,time_interp');
    
    Fs = 1/mean(diff(time_interp));
    NFFT = 2^nextpow2(length(time_interp)); % Next power of 2 from length of y
    Y = fft(data_interp-mean(data_interp),NFFT)/length(time_interp);
    freq = Fs/2*linspace(0,1,NFFT/2+1);
    amp = 2*abs(Y(1:NFFT/2+1));
end