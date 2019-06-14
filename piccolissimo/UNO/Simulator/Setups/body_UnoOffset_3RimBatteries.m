disp('Using body_UnoOffset_3RimBatteries');

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_drag = Cld.Cl_segmented;
Cd_drag = Cld.Cd_segmented;

%% Stator properties
h_s = .02623 - h_cg;
m_s = .172; %stator mass in kg
Is_s = [1900000,0,0;0,1900000,0;0,0,3510000]*1e-9; %Stator inertia at stator cg in flyer frame % rim batteries

%% Drag plate properties
span_offset = .04/2; % span_offset_diameter/2
R_d = span_offset + .165; %Single plate radius (meters) (span_offset_diameter/2 + span)
span_d_nonuniform = span_offset + [-span_offset*1000 0 10 155 160 165]/1000; % UNO V4Base
beta_d_nonuniform = deg2rad(180-[30 30 30 30 0 0]); % UNO V4Base
chord_d_nonuniform = [0 0 80 80 80 80]/1000; % UNO V4Base
B_d = 6; %number of dragplates