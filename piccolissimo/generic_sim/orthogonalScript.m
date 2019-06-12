global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter
% h_r_vals = [-.136 -.068 0 .068 .136];
% h_r_vals = [-.068 0 .068];
% h_r_vals = [.068];
h_r_vals = [0];
% Ir_r_x_vals = [203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2];
% Ir_r_x_vals = [203536*1e-9/2 203536*1e-9 203536*1e-9*2];
Ir_r_x_vals = [203536*1e-9];
% Ir_r_y_vals = [26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2];
% Ir_r_y_vals = [26959*1e-9/2 26959*1e-9 26959*1e-9*2];
Ir_r_y_vals = [26959*1e-9];
Ir_r_z_vals = [182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2];
% Ir_r_z_vals = [182477*1e-9];
% h_s_vals = [-.026*2 -.026 0 .026 .026*2];
% h_s_vals = [-.026 0 .026];
h_s_vals = [0];
% h_s_vals = [-.026];
Is_s_xy_vals = [600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2];
% Is_s_xy_vals = [600000*1e-9/2 600000*1e-9 600000*1e-9*2];
% Is_s_xy_vals = [600000*1e-9];
Is_s_z_vals = [65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2];
% Is_s_z_vals = [65063*1e-9];
h_p_vals = [-.017*2 -.017 0 .017 .017*2];
% h_p_vals = [-.017 0 .017];
% h_p_vals = [0];
% h_p_vals = [.017];
h_d_vals = [-.01*2 -.01 0 .01 .01*2]; %not default spacing because default is so freaking small
% h_d_vals = [-.01 0 .01]; %not default spacing because default is so freaking small
% h_d_vals = [0];
% h_d_vals = [-.001];
% h_d_vals = [.01];
chord_d_vals = [.5 1/1.5 1 1.5 2];
% chord_d_vals = [1];

%MASS?????
%beta?????

setup_flyer_ideal_prop_on_bottom;

%Pick one and use corresponding below
% load F6L5;
load F6L5extra;
% load F1L2F11L5;

% trials = zeros(size(F1L2F11L5,1),19);
% trials = zeros(size(F6L5,1),19);
trials = zeros(size(F6L5extra,1),19);

%h_r_averages = zeros(numParams,numSavedMetrics,paramTrialNum); %remembers each of the
%results for each setting
h_r_averages = zeros(size(h_r_vals,2),9);
%h_r_trials = zeros(numParams,1); %remembers how many times each setting
%has been used
h_r_trials = zeros(size(h_r_vals,2),1);

Ir_r_x_averages = zeros(size(Ir_r_x_vals,2),9);
Ir_r_x_trials = zeros(size(Ir_r_x_vals,2),1);

Ir_r_y_averages = zeros(size(Ir_r_y_vals,2),9);
Ir_r_y_trials = zeros(size(Ir_r_y_vals,2),1);

Ir_r_z_averages = zeros(size(Ir_r_z_vals,2),9);
Ir_r_z_trials = zeros(size(Ir_r_z_vals,2),1);

h_s_averages = zeros(size(h_s_vals,2),9);
h_s_trials = zeros(size(h_s_vals,2),1);

Is_s_xy_averages = zeros(size(Is_s_xy_vals,2),9);
Is_s_xy_trials = zeros(size(Is_s_xy_vals,2),1);

Is_s_z_averages = zeros(size(Is_s_z_vals,2),9);
Is_s_z_trials = zeros(size(Is_s_z_vals,2),1);

h_p_averages = zeros(size(h_p_vals,2),9);
h_p_trials = zeros(size(h_p_vals,2),1);

h_d_averages = zeros(size(h_d_vals,2),9);
h_d_trials = zeros(size(h_d_vals,2),1);

chord_d_averages = zeros(size(chord_d_vals,2),9);
chord_d_trials = zeros(size(chord_d_vals,2),1);

% for trial = 1:size(F1L2F11L5,1)
% for trial = 1:size(F6L5,1)
for trial = 1:size(F6L5extra,1)
% for trial = 35
    trial
    %Figure out which values to test
    %Pick one and use corresponding above
    %Set 1
%         %F1L2F11L5(trial,1) %only has 2 settings
%         h_r_num = F1L2F11L5(trial,2);
%         Ir_r_x_num = F1L2F11L5(trial,3);
%         Ir_r_y_num = F1L2F11L5(trial,4);
%         Ir_r_z_num = F1L2F11L5(trial,5);
%         h_s_num = F1L2F11L5(trial,6);
%         Is_s_xy_num = F1L2F11L5(trial,7);
%         Is_s_z_num = F1L2F11L5(trial,8); 
%         h_p_num = F1L2F11L5(trial,9);
%         h_d_num = F1L2F11L5(trial,10);
%         chord_d_num = F1L2F11L5(trial,11);
%         %F1L2F11L5(trial,12) %has 5 settings
    %Set 2
%         h_r_num = 1;
%         Ir_r_x_num = 1;
%         Ir_r_y_num = 1;
%         Ir_r_z_num = F6L5(trial,1);
%         h_s_num = 1;
%         Is_s_xy_num = F6L5(trial,2);
%         Is_s_z_num = F6L5(trial,3); 
%         h_p_num = F6L5(trial,4);
%         h_d_num = F6L5(trial,5);
%         chord_d_num = F6L5(trial,6);
    %Set 3
        h_r_num = 1;
        Ir_r_x_num = 1;
        Ir_r_y_num = 1;
        Ir_r_z_num = F6L5extra(trial,1);
        h_s_num = 1;
        Is_s_xy_num = F6L5extra(trial,2);
        Is_s_z_num = F6L5extra(trial,3); 
        h_p_num = F6L5extra(trial,4);
        h_d_num = F6L5extra(trial,5);
        chord_d_num = F6L5extra(trial,6);
    
    %Assign values
    h_r = h_r_vals(h_r_num);
    Ir_r_x = Ir_r_x_vals(Ir_r_x_num);
    Ir_r_y = Ir_r_y_vals(Ir_r_y_num);
    Ir_r_z = Ir_r_z_vals(Ir_r_z_num);
    h_s = h_s_vals(h_s_num);
    Is_s_xy = Is_s_xy_vals(Is_s_xy_num);
    Is_s_z = Is_s_z_vals(Is_s_z_num);
    h_p = h_p_vals(h_p_num);
    h_d = h_d_vals(h_d_num);
    chord_d_pre = chord_d_vals(chord_d_num);
    
    chord_d = [.152 .152 .152 .152 .152 .152 .152 .152 .152]*chord_d_pre;

    Ir_r(1,1) = Ir_r_x;
    Ir_r(2,2) = Ir_r_y;
    Ir_r(3,3) = Ir_r_z;
    Is_s(1,1) = Is_s_xy;
    Is_s(2,2) = Is_s_xy;
    Is_s(3,3) = Is_s_z;

    %Calculated Parameters
%     S_s_f = [0 0 h_s];
%     S_r_f = [0 0 h_r];
%     area_p = R_p*R_p*pi;
%     area_d = R_d*R_d*pi;
%     bladeProperties_p = 0;
%     bladeProperties_p(1) = -h_p;
%     bladeProperties_p(2) = R_p;
%     bladeProperties_p(3) = size(beta_p,2);
%     bladeProperties_p = [bladeProperties_p beta_p];
%     bladeProperties_p = [bladeProperties_p chord_p];
%     bladeProperties_d1 = 0;
%     bladeProperties_d1(1) = -h_d;
%     bladeProperties_d1(2) = R_d1;
%     bladeProperties_d1(3) = size(beta_d,2);
%     bladeProperties_d1 = [bladeProperties_d1 beta_d];
%     bladeProperties_d1 = [bladeProperties_d1 chord_d];
%     bladeProperties_d2 = 0;
%     bladeProperties_d2(1) = -h_d;
%     bladeProperties_d2(2) = R_d2;
%     bladeProperties_d2(3) = size(beta_d,2);
%     bladeProperties_d2 = [bladeProperties_d2 beta_d];
%     bladeProperties_d2 = [bladeProperties_d2 chord_d];
%     m = m_s+m_r;
%     I_tot = (m_s*(S_s_f')*S_s_f+m_r*(S_r_f')*S_r_f);

    %Run the simulation
    X = [3.02; 0; 0; 0; 0; 0; 0; -340; 49; .1; 0; 0; 0; 0; 0; 0; 0; -2.524];
    [tout Xout] = yimFlyerLite(X,10);

    %Analyze
    trials(trial,1) = range(Xout(:,10));
    trials(trial,2) = range(Xout(:,11));
    trials(trial,3) = range(sqrt(Xout(:,10).^2. + Xout(:,11).^2));
    trials(trial,4) = std(Xout(:,10));
    trials(trial,5) = std(Xout(:,11));
    trials(trial,6) = std(sqrt(Xout(:,10).^2. + Xout(:,11).^2));
    [peakLoc1, peakMag] = peakfinder(Xout(:,10),trials(trial,1)/6,0,1);
    p1 = polyfit(peakLoc1,peakMag,2);
    f1test = polyval(p1,peakLoc1);
    f1 = polyval(p1,1:size(Xout(:,10)));
%     [peakLoc, peakMag] = peakfinder(Xout(:,11),trials(trial,1)/12,0,1);
%     p2 = polyfit(peakLoc,peakMag,2);
%     f2 = polyval(p2,1:size(Xout(:,10)));
%     [peakLoc, peakMag] = peakfinder(sqrt(Xout(:,10).^2. + Xout(:,11).^2),trials(trial,1)/12,0,1);
%     p3 = polyfit(peakLoc,peakMag,2);
%     f3 = polyval(p3,1:size(Xout(:,10)));
    figure(2);
    hold on;
    plot(tout, f1,'y');
%     plot(tout, f2,'m');
%     plot(tout, f3,'c');
    plot(tout(peakLoc1), f1test,'o');
    pause(.25);
    p1root = p1;
    p1root(3) = p1root(3)-.0368;
    if isreal(min(max(0,roots(p1root)))) && min(max(0,roots(p1root))) >= .5 && size(p1,2)> 2
        trials(trial,7) = tout(round(min(max(0,roots(p1root)))))/100; %/100 to plot nicely!!!!!
    else
        trials(trial,7) = .1;
    end
%     trials(trial,8) = interp1q(tout, f2', .0368);
%     trials(trial,9) = interp1q(tout, f3', .0368);

    %populate
    trials(trial,10) = h_r;
%     h_r_averages(h_r_num,1:9) = h_r_averages(h_r_num,1:9) + trials(trial,1:9);
    h_r_averages(h_r_num,1:9,h_r_trials(h_r_num) + 1) = trials(trial,1:9);
    h_r_trials(h_r_num) = h_r_trials(h_r_num) + 1;

    trials(trial,11) = Ir_r_x;
%     Ir_r_x_averages(Ir_r_x_num,1:9) = Ir_r_x_averages(Ir_r_x_num,1:9) + trials(trial,1:9);
    Ir_r_x_averages(Ir_r_x_num,1:9,Ir_r_x_trials(Ir_r_x_num) + 1) = trials(trial,1:9);
    Ir_r_x_trials(Ir_r_x_num) = Ir_r_x_trials(Ir_r_x_num) + 1;

    trials(trial,12) = Ir_r_y;
%     Ir_r_y_averages(Ir_r_y_num,1:9) = Ir_r_y_averages(Ir_r_y_num,1:9) + trials(trial,1:9);
    Ir_r_y_averages(Ir_r_y_num,1:9,Ir_r_y_trials(Ir_r_y_num) + 1) = trials(trial,1:9);
    Ir_r_y_trials(Ir_r_y_num) = Ir_r_y_trials(Ir_r_y_num) + 1;

    trials(trial,13) = Ir_r_z;
%     Ir_r_z_averages(Ir_r_z_num,1:9) = Ir_r_z_averages(Ir_r_z_num,1:9) + trials(trial,1:9);
    Ir_r_z_averages(Ir_r_z_num,1:9,Ir_r_z_trials(Ir_r_z_num) + 1) = trials(trial,1:9);
    Ir_r_z_trials(Ir_r_z_num) = Ir_r_z_trials(Ir_r_z_num) + 1;

    trials(trial,14) = h_s;
%     h_s_averages(h_s_num,1:9) = h_s_averages(h_s_num,1:9) + trials(trial,1:9);
    h_s_averages(h_s_num,1:9,h_s_trials(h_s_num) + 1) = trials(trial,1:9);
    h_s_trials(h_s_num) = h_s_trials(h_s_num) + 1;

    trials(trial,15) = Is_s_xy;
%     Is_s_xy_averages(Is_s_xy_num,1:9) = Is_s_xy_averages(Is_s_xy_num,1:9) + trials(trial,1:9);
    Is_s_xy_averages(Is_s_xy_num,1:9,Is_s_xy_trials(Is_s_xy_num) + 1) = trials(trial,1:9);
    Is_s_xy_trials(Is_s_xy_num) = Is_s_xy_trials(Is_s_xy_num) + 1;

    trials(trial,16) = Is_s_z;
%     Is_s_z_averages(Is_s_z_num,1:9) = Is_s_z_averages(Is_s_z_num,1:9) + trials(trial,1:9);
    Is_s_z_averages(Is_s_z_num,1:9,Is_s_z_trials(Is_s_z_num) + 1) = trials(trial,1:9);
    Is_s_z_trials(Is_s_z_num) = Is_s_z_trials(Is_s_z_num) + 1;

    trials(trial,17) = h_p;
%     h_p_averages(h_p_num,1:9) = h_p_averages(h_p_num,1:9) + trials(trial,1:9);
    h_p_averages(h_p_num,1:9,h_p_trials(h_p_num) + 1) = trials(trial,1:9);
    h_p_trials(h_p_num) = h_p_trials(h_p_num) + 1;

    trials(trial,18) = h_d;
%     h_d_averages(h_d_num,1:9) = h_d_averages(h_d_num,1:9) + trials(trial,1:9);
    h_d_averages(h_d_num,1:9,h_d_trials(h_d_num) + 1) = trials(trial,1:9);
    h_d_trials(h_d_num) = h_d_trials(h_d_num) + 1;

    trials(trial,19) = chord_d_pre;
%     chord_d_averages(chord_d_num,1:9) = chord_d_averages(chord_d_num,1:9) + trials(trial,1:9);
    chord_d_averages(chord_d_num,1:9,chord_d_trials(chord_d_num) + 1) = trials(trial,1:9);
    chord_d_trials(chord_d_num) = chord_d_trials(chord_d_num) + 1;
end

%process
h_r_means = mean(h_r_averages,3);
h_r_stds = std(h_r_averages,0,3);
% h_r_averages(:,1:9) = h_r_averages(:,1:9)./repmat(h_r_trials,1,9);

Ir_r_x_means = mean(Ir_r_x_averages,3);
Ir_r_x_stds = std(Ir_r_x_averages,0,3);
% Ir_r_x_averages(:,1:9) = Ir_r_x_averages(:,1:9)./repmat(Ir_r_x_trials,1,9);

Ir_r_y_means = mean(Ir_r_y_averages,3);
Ir_r_y_stds = std(Ir_r_y_averages,0,3);
% Ir_r_y_averages(:,1:9) = Ir_r_y_averages(:,1:9)./repmat(Ir_r_y_trials,1,9);

Ir_r_z_means = mean(Ir_r_z_averages,3);
Ir_r_z_stds = std(Ir_r_z_averages,0,3);
% Ir_r_z_averages(:,1:9) = Ir_r_z_averages(:,1:9)./repmat(Ir_r_z_trials,1,9);

h_s_means = mean(h_s_averages,3);
h_s_stds = std(h_s_averages,0,3);
% h_s_averages(:,1:9) = h_s_averages(:,1:9)./repmat(h_s_trials,1,9);

Is_s_xy_means = mean(Is_s_xy_averages,3);
Is_s_xy_stds = std(Is_s_xy_averages,0,3);
% Is_s_xy_averages(:,1:9) = Is_s_xy_averages(:,1:9)./repmat(Is_s_xy_trials,1,9);

Is_s_z_means = mean(Is_s_z_averages,3);
Is_s_z_stds = std(Is_s_z_averages,0,3);
% Is_s_z_averages(:,1:9) = Is_s_z_averages(:,1:9)./repmat(Is_s_z_trials,1,9);

h_p_means = mean(h_p_averages,3);
h_p_stds = std(h_p_averages,0,3);
% h_p_averages(:,1:9) = h_p_averages(:,1:9)./repmat(h_p_trials,1,9);

h_d_means = mean(h_d_averages,3);
h_d_stds = std(h_d_averages,0,3);
% h_d_averages(:,1:9) = h_d_averages(:,1:9)./repmat(h_d_trials,1,9);

chord_d_means = mean(chord_d_averages,3);
chord_d_stds = std(chord_d_averages,0,3);
% chord_d_averages(:,1:9) = chord_d_averages(:,1:9)./repmat(chord_d_trials,1,9);

%generate colors
ColorSet = varycolor(9);
n = 1:7;

%plot
figure(9);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(h_r_vals,h_r_means,'LineWidth', 2);
plot(h_r_vals,h_r_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('h_r');

figure(10);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(Ir_r_x_vals,Ir_r_x_means,'LineWidth', 2);
plot(Ir_r_x_vals,Ir_r_x_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('Irr_x');

figure(11);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(Ir_r_y_vals,Ir_r_y_means,'LineWidth', 2);
plot(Ir_r_y_vals,Ir_r_y_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('Irr_y');

figure(12);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(Ir_r_z_vals,Ir_r_z_means,'LineWidth', 2);
plot(Ir_r_z_vals,Ir_r_z_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('Irr_z');

figure(13);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(h_s_vals,h_s_means,'LineWidth', 2);
plot(h_s_vals,h_s_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('h_s');

figure(14);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(Is_s_xy_vals,Is_s_xy_means,'LineWidth', 2);
plot(Is_s_xy_vals,Is_s_xy_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('Iss_x_y');

figure(15);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(Is_s_z_vals,Is_s_z_means,'LineWidth', 2);
plot(Is_s_z_vals,Is_s_z_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('Iss_z');

figure(16);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(h_p_vals,h_p_means,'LineWidth', 2);
plot(h_p_vals,h_p_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('h_p');

figure(17);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(h_d_vals,h_d_means,'LineWidth', 2);
plot(h_d_vals,h_d_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('h_d');

figure(18);
cla; hold on; grid on; 
set(gca, 'ColorOrder', ColorSet);
plot(chord_d_vals,chord_d_means,'LineWidth', 2);
plot(chord_d_vals,chord_d_stds,'LineWidth', 2,'LineStyle',':');
legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
title('chord_d');
