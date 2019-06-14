disp('Using body_UnoV2_3RimBatteries');

%% Airfoil data
Cld = load('NACA0012_RE360000_segmented_HR.mat');
Cl_drag = Cld.Cl_segmented;
Cd_drag = Cld.Cd_segmented;

%% Stator properties
h_s = .03535 - h_cg;
m_s = .177; %stator mass in kg
Is_s = [1630000,0,0;0,1630000,0;0,0,3020000]*1e-9; %Stator inertia at stator cg in flyer frame % rim batteries

%% Drag plate properties
span_offset = .04/2; % span_offset_diameter/2
R_d = span_offset + .155; %Single plate radius (meters)
span_d_nonuniform = [0 20 51 128.75 175]/1000;
beta_d_nonuniform = deg2rad(180-[90 70 53 (27.5) (21)]); % UNO V2
chord_d_nonuniform = [0 .035 .052 .085 .117]; % UNO V2
B_d = 6; %number of dragplates