disp('Using body_UnoV3_3RimBatteries');

%% Airfoil data
Cld = load('NACA009_RE40000_fab_segmented_HR.mat');
Cl_drag = Cld.Cl_segmented;
Cd_drag = Cld.Cd_segmented;

%% Stator properties
h_s = .02963 - h_cg;
m_s = .177;%154.72/1000; %stator mass in kg
Is_s = [1840000,0,0;0,1840000,0;0,0,3450000]*1e-9; %Stator inertia at stator cg in flyer frame % rim batteries

%% Drag plate properties
span_offset = .04/2; % span_offset_diameter/2
R_d = span_offset + .165; %Single plate radius (meters) (span_offset_diameter/2 + span)
span_d_nonuniform = span_offset + [-span_offset*1000 0 33 115.5 165]/1000;
beta_d_nonuniform = deg2rad(180-[90 70 53 27.5 21]); % UNO V3
chord_d_nonuniform = [0 35 52 95 120]/1000; % UNO V3
B_d = 6; %number of dragplates