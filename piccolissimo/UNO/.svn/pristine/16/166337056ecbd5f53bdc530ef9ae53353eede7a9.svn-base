disp('Using propeller_gemfan_12x4_7LL_snipped');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -322.1825; 46.6698; 0; 0; 0; 0; 0; 0; 0; 0; -1.411];

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

rise = [0 0 1.34 5.66 8.99 8.69 8.77 7.62 5.16 4.37]/1000;
run = [0 0 12.04 15.81 25.61 29.76 31.2 28.99 22.16 19.03]/1000;
span = .025 + [0 8 10 20 40 60 80 100 120 127]/1000;

m_r = .039; %TODO::find real number, I'm an estimate %rotor mass in kg
Ir_r = [ 8550,0,0;0,8300+66500,0;0,0, 8300+66500]*1e-9; %TODO::find real number, I'm an estimate %Rotor inertia at rotor cg in flyer frame 
R_p = max(span); %Single blade radius (meters) %12x4.7
beta_p_nonuniform = atan(rise./run);
beta_p_nonuniform(isnan(beta_p_nonuniform)) = 0;
chord_p_nonuniform = sqrt(rise.^2+run.^2);
span_p_nonuniform = span;
B_p = 2; %number of propeller blades