disp('Using propeller_LLFlap_10_deg_0_twist_Quanum2206');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -240.5759; 39.5377; 0; 0; 0; 0; 0; 0; 0; 0; -1.815];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = .02648;
Ir_r = [ 3110,0,0;0,85000,0;0,0,85000]*1e-9; % no lg %Rotor inertia at rotor cg in flyer frame 
R_p = .157; %Single blade radius (meters) %
beta_p_nonuniform = deg2rad(10)*[1 1 1 1 1 1];
chord_p_nonuniform = .0195*[1 1 1 1 1 1];
span_p_nonuniform = linspace(0,R_p,length(beta_p_nonuniform));
B_p = 2; %number of propeller blades