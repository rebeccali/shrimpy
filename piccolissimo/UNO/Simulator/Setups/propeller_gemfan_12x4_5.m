disp('Using propeller_gemfan_12x4.5');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -322.1825; 46.6698; 0; 0; 0; 0; 0; 0; 0; 0; -1.411];

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = 10.93/1000 + 12.09/1000; %=.02302; %propeller mass + dys quanum 2206 rotor mass in kg
Ir_r = [ (376+885.08),0,0;0,(65500+885.08),0;0,0, (65500+1336.33)]*1e-9; %propeller inertia + dys quanum 2206 rotor %Rotor inertia at rotor cg in flyer frame

rise = [1.34 5.66 8.99 8.69 8.77 7.62 5.16 4.37 3.42 1.31]/1000; %in meters
run = [12.04 15.81 25.61 29.76 31.2 28.99 22.16 19.03 14.4 10.44]/1000;
span = [10 20 40 60 80 100 120 127 135 138.5]/1000;

R_p = .150; %Single blade radius (meters) %12x4.5
beta_p_nonuniform = atan([0 9.25/13.2 8/22.75 7.4/26.5 5.7/25.5 3.5/20.5 2/11]); %12x4.5
chord_p_nonuniform = [.0127 .0161 .0241 .0275 .0261 .0208 .0112]; %12x4.5
span_p_nonuniform = linspace(0,R_p,length(beta_p_nonuniform));
B_p = 2; %number of propeller blades