disp('Using UnoControl');
global kpp kdp kip r_b kff_pitch kpv kdv kiv kp_pitch kd_pitch kdd_pitch pitch_clamp kpomg kdomg kiomg v_clamp psi_0 waypoints target_size flyer_ctrl_fcn

%% Battery setup
v_clamp = 11.1;
r_b = 0; % battery resistance in ohms

%% Desired waypoints
%     waypoints = [1 0 0; 1 0 1; 1 1 1; 0 1 1; 0 1 0; 0 0 0];
%     waypoints = [0 0 0; 0 0 1; 1 0 1; 1 0 0; 1 1 0; 1 1 1; 0 1 1; 0 1 0; 0 0 0];
%     waypoints = [10 0 0; 20 0 0];
    waypoints = [0 0 0];
%     waypoints = [999 0 0];
    target_size = .15;

%% Gains
    % Voltage gains
%     kpv = 1;
%     kdv = 0;%K_t_m*.01;
%     kiv = 0;%K_t_m;
    kpp = K_t_m/3.7;
    kdp = 0;
    kip = 0;

    % Speed gains
    kpomg = 5000;
    kdomg = -.1*kpomg;
    kiomg = 0;%kpomg*.05;
    
    % Pulsing gains
    kp_pitch = 5; % degrees/m
    kd_pitch = 1; % degrees/(m/s)
    kdd_pitch = 0; % degrees/(body radian)
    kff_pitch = 0; % degrees/(m/s^3)
    pitch_clamp = 10;
    psi_0 = -atan2(3.141,-.6264); % phase shift to make vehicle tilt/move in x direction (tilt about -y)
    
%% Controller
    flyer_ctrl_fcn = @UnoControlLoop;