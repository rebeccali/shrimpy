disp('Using propeller_gemfan_11x4.7');

Xbase = [nan; 0; 0; 0; 0; 0; 0; -381.7049; 45.7068; 0; 0; 0; 0; 0; 0; 0; 0; -1.242];

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_prop = Cld.Cl_segmented;
Cd_prop = Cld.Cd_segmented;

m_r = 9.71/1000 + 12.09/1000; %propeller mass + dys quanum 2206 rotor mass in kg
Ir_r = [(341+885.08),0,0;0,(58200+885.08),0;0,0, (58200+1336.33)]*1e-9;%Rotor inertia at rotor cg in flyer frame 
R_p = .140; %Single blade radius (meters) %11x4.7
beta_p_nonuniform = atan([0 4.24/13.72 7.81/19.08 9.62/24.47 9.54/28.47 8.33/30.83 6.84/31.04 5.3/29.1 3.73/24.8 2.1/18.21 .34/4.89]); %11x4.7
chord_p_nonuniform = [.015 .01437 .02062 .0263 .03003 .03193 .03178 .02958 .02508 .01833 .0049];%11x4.7
span_p_nonuniform = linspace(0,R_p,length(beta_p_nonuniform));
B_p = 2; %number of propeller blades