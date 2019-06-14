disp('Using propeller_gemfan_10x4.5');

% Not finished
Xbase = [nan; 0; 0; 0; 0; 0; 0; -419.127632404831; 45.4149756916852; 0; 0; 0; 0; 0; 0; 0; 0; -1.242];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = 7.73/1000 + 12.09/1000; %=.02302; %rotor mass + dys quanum 2206 rotor mass in kg
Ir_r = [(165+885.08),0,0;0,(36700+885.08),0;0,0, (36700+1336.33)]*1e-9; %propeller inertia + dys quanum 2206 rotor %Rotor inertia at rotor cg in flyer frame 

rise = [7.5 8.5 8.5 8.5 5 3.5 2 2]/1000;
run = [5 12 21 26 26 22 14 8]/1000;
span = [5 20 40 60 80 100 120 125]/1000;

R_p = max(span); %Single blade radius (meters) 
beta_p_nonuniform = atan(rise./run);
beta_p_nonuniform(isnan(beta_p_nonuniform)) = 0;
chord_p_nonuniform = sqrt(rise.^2+run.^2);
span_p_nonuniform = span;
B_p = 2; %number of propeller blades