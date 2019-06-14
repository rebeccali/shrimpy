disp('Using PiccolissimoControl');
global r_b kppsi kdpsi kipsi kpp kdp kip kpomg kdomg kiomg v_clamp psi_0 waypoints flyer_ctrl_fcn


%% Battery setup
    v_clamp = 3.7;
    r_b = 2;

%% Desired waypoints
%     waypoints = [0 0 0];
    waypoints = [999 0 0];
%     waypoints = [999 -800 0];
    target_size = .05;

%% Gains
    % Voltage gains
%     kpv = K_t_m;
%     kdv = 0;%K_t_m*.01;
%     kiv = 0;%K_t_m;
    kpp = K_t_m/3.7;
    kdp = 0;
    kip = 0;

    % Speed gains
    kpomg = -1000;
    kdomg = -.2*kpomg;
    kiomg = 0;%kpomg*.05;
    
    % Pulsing gains
    kppsi = 0;
    kdpsi = 0;
    kipsi = 0;
    psi_0 = pi/4+pi+pi/16+1/32; % phase shift to make vehicle tilt/move in x direction (tilt about -y)
    
%% Controller
    flyer_ctrl_fcn = @PiccolissimoControlLoop;