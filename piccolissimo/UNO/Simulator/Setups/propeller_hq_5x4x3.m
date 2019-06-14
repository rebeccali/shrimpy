disp('Using propeller_hq_5x4x3');

% Not finished
Xbase = [nan; 0; 0; 0; 0; 0; 0; -419.127632404831; 45.4149756916852; 0; 0; 0; 0; 0; 0; 0; 0; -1.242];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = 8.95/1000 + 12.09/1000; %=.02302; %propeller mass + dys quanum 2206 rotor mass in kg
Ir_r = [(210+885.08),0,0;0,(35200+885.08),0;0,0, (35200+1336.33)]*1e-9; %propeller inertia + dys quanum 2206 rotor %Rotor inertia at rotor cg in flyer frame 

rise = [5 6 5 1 0.5]/1000; %in meters
run = [14 14 13 8 4]/1000;
span = [6 20 40 60 76]/1000;

R_p = max(span); %Single blade radius (meters) 
beta_p_nonuniform = atan(rise./run);
beta_p_nonuniform(isnan(beta_p_nonuniform)) = 0;
chord_p_nonuniform = sqrt(rise.^2+run.^2);
span_p_nonuniform = span;
B_p = 2; %number of propeller blades

