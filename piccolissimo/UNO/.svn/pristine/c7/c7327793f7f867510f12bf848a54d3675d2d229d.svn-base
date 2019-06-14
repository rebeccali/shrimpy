disp('Using propeller_6in_14_deg_0_twist_Quanum2206');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -936.4877; 22.4905; 0; 0; 0; 0; 0; 0; 0; 0; -4.862];

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = .02648;
Ir_r = [ 3110,0,0;0,85000,0;0,0,85000]*1e-9; % no lg %Rotor inertia at rotor cg in flyer frame 
R_p = .153/2; %Single blade radius (meters) %
beta_p_nonuniform = deg2rad(14)*[1 1 1 1 1 1];
chord_p_nonuniform = .0195*[1 1 1 1 1 1];
span_p_nonuniform = linspace(0,R_p,length(beta_p_nonuniform));
B_p = 2; %number of propeller blades