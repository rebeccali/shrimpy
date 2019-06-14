disp('Using propeller_unk_8x4.5');

% Not finished
Xbase = [nan; 0; 0; 0; 0; 0; 0; -419.127632404831; 45.4149756916852; 0; 0; 0; 0; 0; 0; 0; 0; -1.242];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = 5.93/1000 + 12.09/1000; %=.02302; %propeller mass + dys quanum 2206 rotor mass in kg
Ir_r = [ (93.3+885.08),0,0;0,(18000+885.08),0;0,0, (18000+1336.33)]*1e-9; %propeller inertia + dys quanum 2206 rotor %Rotor inertia at rotor cg in flyer frame 

rise = [7 8.5 7.5 7 4.5 2.5 2.5]/1000; %in meters
run = [10 12 19 22 20 11 8]/1000;
span = [6 20 40 60 80 100 101.5]/1000;

R_p = max(span); %Single blade radius (meters) 
beta_p_nonuniform = atan(rise./run);
beta_p_nonuniform(isnan(beta_p_nonuniform)) = 0;
chord_p_nonuniform = sqrt(rise.^2+run.^2);
span_p_nonuniform = span;
B_p = 2; %number of propeller blades

