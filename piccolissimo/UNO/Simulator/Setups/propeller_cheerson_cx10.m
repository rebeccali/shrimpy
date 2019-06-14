disp('Using propeller_cheerson_cx10');

prop_adj = 0;%-deg2rad(4.5)

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

R_p = .0146; %Cheerson %Single blade radius (meters) 
span_p_nonuniform = 0:.00292:R_p;
beta_p_nonuniform = prop_adj-atan([0 .56/3.17 .72/3.91 .75/4.06 .63/3.48 .15/.89])'; %Cheerson taken every 2.92mm
chord_p_nonuniform = [0 3.22/1000 3.98/1000 4.12/1000 3.54/1000 .90/1000]'; %Cheerson taken every 2.92mm

m_r = .25/1000; %rotor mass in kg
Ir_r = [3.95, 0, 0; 0, 6.31, 0; 0, 0, 3.21]*1e-9; %Rotor inertia at rotor cg in flyer frame 
B_p = 2; %number of propeller blades