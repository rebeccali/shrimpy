disp('Using propeller_fixed_14_deg_0_twist');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -262.3; 42.55; 0; 0; 0; 0; 0; 0; 0; 0; -4.862];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = .046;
Ir_r = [ 8550,0,0;0,8300+71000,0;0,0, 8300+71300]*1e-9; % no lg %Rotor inertia at rotor cg in flyer frame 
R_p = .157; %Single blade radius (meters) %
beta_p_nonuniform = deg2rad(14)*[1 1 1 1 1 1];
chord_p_nonuniform = .0195*[1 1 1 1 1 1];
span_p_nonuniform = linspace(0,R_p,length(beta_p_nonuniform));
B_p = 2; %number of propeller blades