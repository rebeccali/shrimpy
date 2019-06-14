disp('Using body_Piccolissimo_V11');

drag_adj = -deg2rad(2);
inertia_adj = 2.5;

%% Airfoil data
Cld = load('FlatPlateLentink_RE14000_segmented_HR.mat');
Cl_drag = Cld.Cl_segmented;
Cd_drag = Cld.Cd_segmented;

%% Stator properties
h_s = -2.69/1000 - h_cg;
m_s = (4.47-.25)/1000;   %stator mass in kg
Is_s = inertia_adj*[383, 0, 0; 0, 439, 0; 0, 0, 697]*1e-9; %Stator inertia at stator cg in flyer frame % rim batteries

%% Drag plate properties
span_d_nonuniform = [0 3.18 5.31 9.58 15.91 19.16]/1000;
rise_d_nonuniform = [13.5 13.43 12.44 10.63 7.34 3.95]/1000;
run_d_nonuniform = [0 1.39 3.06 5.80 5.86 3.63]/1000;
beta_d_nonuniform = drag_adj + atan(rise_d_nonuniform./run_d_nonuniform);
beta_d_nonuniform(isnan(beta_d_nonuniform)) = 0;
chord_d_nonuniform = sqrt(rise_d_nonuniform.^2+run_d_nonuniform.^2);
R_d = max(span_d_nonuniform); %Single plate radius (meters)
B_d = 6; %number of dragplates