global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d1 R_d2 beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter



% valsMat = [-.136 -.068 0 .068 .136; ... %h_r_vals = [-.136 -.068 0 .068 .136];
%     203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2; ... %Ir_r_x_vals = [203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2];
%     26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2; ... % Ir_r_y_vals = [26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2];
%     182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2; ... % Ir_r_z_vals = [182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2];
%     -.026*2 -.026 0 .026 .026*2; ... % h_s_vals = [-.026*2 -.026 0 .026 .026*2];
%     600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2; ... %Is_s_xy_vals = [600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2];
%     65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2; ... %Is_s_z_vals = [65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2];
%     -.017*2 -.017 0 .017 .017*2; ... %h_p_vals = [-.017*2 -.017 0 .017 .017*2];
%     -.01*2 -.01 0 .01 .01*2; ... %h_d_vals = [-.01*2 -.01 0 .01 .01*2]; %not default spacing because default is so freaking small
%     .5 1/1.5 1 1.5 2; ... %chord_d_vals = [.5 1/1.5 1 1.5 2];
%     0 .1 nan nan nan; ... %beta
%     ];

% valsMat = [0 nan nan; ... %h_r_vals = [-.136 -.068 0 .068 .136];
%     203536*1e-9 nan nan; ... %Ir_r_x_vals = [203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2];
%     26959*1e-9 nan nan; ... % Ir_r_y_vals = [26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2];
%     182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5; ... % Ir_r_z_vals = [182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2];
%     0 nan nan; ... % h_s_vals = [-.026*2 -.026 0 .026 .026*2];
%     600000*1e-9 nan nan; ... %Is_s_xy_vals = [600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2];
%     65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5; ... %Is_s_z_vals = [65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2];
%     -.017 0 .017; ... %h_p_vals = [-.017*2 -.017 0 .017 .017*2];
%     -.01 0 .01; ... %h_d_vals = [-.01*2 -.01 0 .01 .01*2]; %not default spacing because default is so freaking small
%     1/1.5 1 1.5; ... %chord_d_vals = [.5 1/1.5 1 1.5 2];
%     -.15 0 .15; ... %beta
%     ];

%4
% valsMat = [0 nan nan nan nan nan nan nan nan; ... %h_r_vals = [-.136 -.068 0 .068 .136];
%     203536*1e-9 nan nan nan nan nan nan nan nan; ... %Ir_r_x_vals = [203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2];
%     26959*1e-9 nan nan nan nan nan nan nan nan; ... % Ir_r_y_vals = [26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2];
%     182477*1e-9/1.5 182477*1e-9*1.5 nan nan nan nan nan nan nan; ... % Ir_r_z_vals = [182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2];
%     0 nan nan nan nan nan nan nan nan; ... % h_s_vals = [-.026*2 -.026 0 .026 .026*2];
%     600000*1e-9 nan nan nan nan nan nan nan nan; ... %Is_s_xy_vals = [600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2];
%     65063*1e-9/1.5 65063*1e-9*1.5 nan nan nan nan nan nan nan; ... %Is_s_z_vals = [65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2];
%     -.017 .017 nan nan nan nan nan nan nan; ... %h_p_vals = [-.017*2 -.017 0 .017 .017*2];
%     -.01 .01 nan nan nan nan nan nan nan; ... %h_d_vals = [-.01*2 -.01 0 .01 .01*2]; %not default spacing because default is so freaking small
%     1/2 1/1.5 1 1.5 2 nan nan nan nan; ... %chord_d_vals = [.5 1/1.5 1 1.5 2];
%     -.2 -.15 -.1 -.05 0 .05 .1 .15 .2; ... %beta
%     ];

%V1
% valsMat = [0 nan; ...
%     203536*1e-9 nan; ...
%     26959*1e-9 nan; ...
%     182477*1e-9/1.5 182477*1e-9*1.5; ...
%     0 nan; ...
%     600000*1e-9 nan; ...
%     65063*1e-9/1.5 65063*1e-9*1.5; ...
%     -.017 .017; ...
%     -.01 .01; ...
%     1/1.5 1.5; ...
%     -.05 .05; ...
%     ];

%V1.5 (corrected V1)
% valsMat = [.068 nan; ... %68mm
%     26959*1e-9 nan; ...
%     203536*1e-9 nan; ...
%     182477*1e-9/1.5 182477*1e-9*1.5; ...
%     -.026 nan; ... %-26mm
%     600000*1e-9 nan; ...
%     66438*1e-9/1.5 66438*1e-9*1.5; ...
%     -.085 .085; ... %84.76mm
%     -.026 .026; ...%-26mm
%     1/1.5 1.5; ...
%     -.05 .05; ...
%     ];

%V2
% valsMat = [.04823 nan; ...%h_r_vals
%     23041.57*1e-9 nan; ...%Ir_r_x_vals
%     199660.67*1e-9 nan; ...%Ir_r_y_vals
%     182456.13*1e-9/1.5 182456.13*1e-9*1.5; ...%Ir_r_z_vals
%     -.03753 nan; ...%h_s_vals
%     500000*1e-9 nan; ...%Is_s_xy_vals
%     106007.57*1e-9/1.5 106007.57*1e-9*1.5; ...%Is_s_z_vals
%     -.0664 .0664; ...%h_p_vals
%     -.08892 .08892; ...%h_d_vals
%     1/1.5 1.5; ...%chord_d_vals
%     -.05 .05; ...%beta
%     ];

%V2.1
valsMat = [.04823 nan nan nan; ...%h_r_vals
    23041.57*1e-9 nan nan nan; ...%Ir_r_x_vals
    199660.67*1e-9 nan nan nan; ...%Ir_r_y_vals
    182456.13*1e-9 nan nan nan; ...%Ir_r_z_vals
    -.03753 nan nan nan; ...%h_s_vals
    500000*1e-9 nan nan nan; ...%Is_s_xy_vals
    106007.57*1e-9 nan nan nan; ...%Is_s_z_vals
    .0664 nan nan nan; ...%h_p_vals
    -.01 -.0051 nan nan; ...%h_d_vals
    1 nan nan nan; ...%chord_d_vals
    0 nan nan nan; ...%beta
    ];

% valsMat = [.04823 nan nan; ...%h_r_vals
%     23041.57*1e-9 nan nan; ...%Ir_r_x_vals
%     199660.67*1e-9 nan nan; ...%Ir_r_y_vals
%     182456.13*1e-9/2 182456.13*1e-9 182456.13*1e-9*2; ...%Ir_r_z_vals
%     -.03753 nan nan; ...%h_s_vals
%     500000*1e-9 nan nan; ...%Is_s_xy_vals
%     106007.57*1e-9/2 106007.57*1e-9 106007.57*1e-9*2; ...%Is_s_z_vals
%     .0664/2 .0664 .0664*1.5; ...%h_p_vals
%     -.08892/2 -.08892 -.08892*1.5; ...%h_d_vals
%     1/1.5 1 1.5; ...%chord_d_vals
%     -.05 0 .05; ...%beta
%     ];
    
%MASS?????

% setup_flyer_ideal_prop_on_bottom;
setup_flyerV2_prop_on_bottom

number_of_metrics = 6;

orderMat = makeMat(sum(~isnan(valsMat),2)');

trials = zeros(size(orderMat,1),number_of_metrics+size(orderMat,2));

trials_ordered = zeros(size(orderMat,2),max(max(orderMat)),number_of_metrics); %num variables, max settings per variable, num metrics, num times this setting is seen (grows dynamically)

trials_num = zeros(size(valsMat));

size(orderMat,1)

for trial = 1:size(orderMat,1)
    trial
   
    %Assign values
    %xxx_num = orderMat(trial,x);
    %xxx = valsMat(xxx_num,x);
    %or
    %xxx = valsMat(orderMat(trial,x),x);
    h_r = valsMat(1,orderMat(trial,1));
    Ir_r_x = valsMat(2,orderMat(trial,2));
    Ir_r_y = valsMat(3,orderMat(trial,3));
    Ir_r_z = valsMat(4,orderMat(trial,4));
    h_s = valsMat(5,orderMat(trial,5));
    Is_s_xy = valsMat(6,orderMat(trial,6));
    Is_s_z = valsMat(7,orderMat(trial,7));
    h_p = valsMat(8,orderMat(trial,8));
    h_d = valsMat(9,orderMat(trial,9));
    chord_d_pre = valsMat(10,orderMat(trial,10));
    beta_d_pre = valsMat(11,orderMat(trial,11));
    
    chord_d = [.152 .152 .152 .152 .152 .152 .152 .152 .152]*chord_d_pre;
    beta_d = [pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2] + beta_d_pre;

    Ir_r(1,1) = Ir_r_x;
    Ir_r(2,2) = Ir_r_y;
    Ir_r(3,3) = Ir_r_z;
    Is_s(1,1) = Is_s_xy;
    Is_s(2,2) = Is_s_xy;
    Is_s(3,3) = Is_s_z;

    %Run the simulation
%     X = [3.02; 0; 0; 0; 0; 0; 0; -340; 49; .1; 0; 0; 0; 0; 0; 0; 0; -2.524];
    X = [2.70; 0; 0; 0; 0; 0; 0; -285.8; 20.53; .1; 0; 0; 0; 0; 0; 0; 0; -1.70];
    [tout Xout] = yimFlyerLite(X,10,'');
    
    sz = size(tout);
    orientation = zeros(sz(1),3);
    for i = 1:sz(1)
        angles(1) = Xout(i,10);
        angles(2) = Xout(i,11);
        angles(3) = Xout(i,12);
        Leb = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        orientation(i,:) = (Leb*[0 0 1]')';
    end

    %Analyze
    trials(trial,1) = range(Xout(:,10));
    trials(trial,2) = range(Xout(:,11));
    trials(trial,3) = range(sqrt(Xout(:,10).^2. + Xout(:,11).^2));
    trials(trial,4) = 1-min(orientation(:,3));
    trials(trial,5) = trapz(abs(1-orientation(:,3)))./100;
%     trials(trial,5) = std(Xout(:,11));
%     trials(trial,6) = std(sqrt(Xout(:,10).^2. + Xout(:,11).^2));
    filterAmount = 100;
    filterOut = filter(ones(1,filterAmount)/filterAmount,1,orientation(:,3));
    [peakLoc1, peakMag] = peakfinder(filterOut(filterAmount:end),trials(trial,4)/25,1,-1);
    if(size(peakLoc1,1) < 3)
        [peakLoc1, peakMag] = peakfinder(filterOut(filterAmount:end),trials(trial,4)/50,1,-1);
    end
    if size(peakLoc1,1) > 2
        p1 = polyfit(peakLoc1,peakMag,2);
        f1test = polyval(p1,peakLoc1);
        f1 = polyval(p1,1:size(orientation(:,3)));

        figure(6);
        hold on;
        plot(tout, f1,'y');
        plot(tout,filterOut,'m');
        plot(tout(peakLoc1), f1test,'o');

        p1root = p1;
    %     p1root(3) = p1root(3)-.0368;
        p1root(3) = p1root(3)-.9982;
        fitRoots = min(min(max(0,roots(p1root))),length(tout));

        if size(p1,2)> 2 && ~isempty(fitRoots) && isreal(fitRoots) &&  fitRoots >= .5
            trials(trial,6) = tout(round(fitRoots))/100; %/100 to plot nicely!!!!!
            plot(tout(round(fitRoots)),.9982,'x');
        else
            trials(trial,6) = .1;
            plot(10,.9982,'x');
        end
    else %TODO::Fix me, I think I return 0 always!!!
%         underTime = orientation(:,3) > .9982;
%         crossOvers = diff(underTime);
%         lastBlipIndex = find(crossOvers ~= 0, 1, 'last' );
%         if crossOvers(lastBlipIndex) > 0
%             trials(trial,6) = tout(lastBlipIndex)/100;
%             plot(tout(lastBlipIndex),.9982,'x');
%         else
            trials(trial,6) = .1;
            plot(10,.9982,'x');
%         end
    end
    axis([0 10 .995 1]);
    pause(.25);
%     trials(trial,8) = interp1q(tout, f2', .0368);
%     trials(trial,9) = interp1q(tout, f3', .0368);

    %populate
    for i = 1:size(orderMat,2)
        trials(trial,number_of_metrics+i) = valsMat(i,orderMat(trial,i));
        trials_ordered(i,orderMat(trial,i),1:number_of_metrics,trials_num(i,orderMat(trial,i))+1) = trials(trial,1:number_of_metrics);
        trials_num(i,orderMat(trial,i)) = trials_num(i,orderMat(trial,i)) + 1;
    end
end

%process
means = zeros(size(orderMat,2),max(max(orderMat)),number_of_metrics);
stds = zeros(size(orderMat,2),max(max(orderMat)),number_of_metrics);

numSettings = sum(~isnan(valsMat),2);
 
for i = 1:size(valsMat,1)
    means(i,1:numSettings(i),:) = mean(trials_ordered(i,1:numSettings(i),:,1:trials_num(i,1)),4); %only works if all settings performed same amount
    stds(i,1:numSettings(i),:) = std(trials_ordered(i,1:numSettings(i),:,1:trials_num(i,1)),0,4); %only works if all settings performed same amount
end

%generate colors
ColorSet = varycolor(number_of_metrics);
titles = {'h_r' 'Irr_x' 'Irr_y' 'Irr_z' 'h_s' 'Iss_x_y' 'Iss_z' 'h_p' 'h_d' 'chord_d' 'beta'};

%plot
for i = 1:size(valsMat,1)
    if numSettings(i) > 1
        figure(i+9);
        cla; hold on; grid on; 
        set(gca, 'ColorOrder', ColorSet);
        plot(valsMat(i,1:numSettings(i)),reshape(means(i,1:numSettings(i),:),numSettings(i),number_of_metrics));
        plot(valsMat(i,1:numSettings(i)),reshape(stds(i,1:numSettings(i),:),numSettings(i),number_of_metrics),'-.');
        title(titles(i));
        legend('phi range','theta range', 'angle range', 'range', 'int/100', 'time const/100');
    end
end

% figure(9);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(h_r_vals,h_r_means,'LineWidth', 2);
% plot(h_r_vals,h_r_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('h_r');
% 
% figure(10);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(Ir_r_x_vals,Ir_r_x_means,'LineWidth', 2);
% plot(Ir_r_x_vals,Ir_r_x_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('Irr_x');
% 
% figure(11);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(Ir_r_y_vals,Ir_r_y_means,'LineWidth', 2);
% plot(Ir_r_y_vals,Ir_r_y_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('Irr_y');
% 
% figure(12);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(Ir_r_z_vals,Ir_r_z_means,'LineWidth', 2);
% plot(Ir_r_z_vals,Ir_r_z_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('Irr_z');
% 
% figure(13);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(h_s_vals,h_s_means,'LineWidth', 2);
% plot(h_s_vals,h_s_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('h_s');
% 
% figure(14);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(Is_s_xy_vals,Is_s_xy_means,'LineWidth', 2);
% plot(Is_s_xy_vals,Is_s_xy_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('Iss_x_y');
% 
% figure(15);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(Is_s_z_vals,Is_s_z_means,'LineWidth', 2);
% plot(Is_s_z_vals,Is_s_z_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('Iss_z');
% 
% figure(16);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(h_p_vals,h_p_means,'LineWidth', 2);
% plot(h_p_vals,h_p_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('h_p');
% 
% figure(17);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(h_d_vals,h_d_means,'LineWidth', 2);
% plot(h_d_vals,h_d_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('h_d');
% 
% figure(18);
% cla; hold on; grid on; 
% set(gca, 'ColorOrder', ColorSet);
% plot(chord_d_vals,chord_d_means,'LineWidth', 2);
% plot(chord_d_vals,chord_d_stds,'LineWidth', 2,'LineStyle',':');
% legend('phi avg range','theta avg range', 'angle avg range', 'phi avg std','theta avg std', 'angle avg std', 'phi avg tau','phi range std','theta range std', 'angle range std', 'phi std std','theta std std', 'angle std std', 'phi tau std');
% title('chord_d');
