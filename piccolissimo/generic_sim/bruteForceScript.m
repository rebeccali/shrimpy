global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter
% h_r_vals = [-.136 -.068 0 .068 .136];
% h_r_vals = [-.068 0 .068];
h_r_vals = [.068];
% h_r_vals = [0];
% Ir_r_x_vals = [203536*1e-9/2 203536*1e-9/1.5 203536*1e-9 203536*1e-9*1.5 203536*1e-9*2];
Ir_r_x_vals = [203536*1e-9/2 203536*1e-9 203536*1e-9*2];
% Ir_r_x_vals = [203536*1e-9];
% Ir_r_y_vals = [26959*1e-9/2 26959*1e-9/1.5 26959*1e-9 26959*1e-9*1.5 26959*1e-9*2];
Ir_r_y_vals = [26959*1e-9/2 26959*1e-9 26959*1e-9*2];
% Ir_r_y_vals = [26959*1e-9];
Ir_r_z_vals = [182477*1e-9/2 182477*1e-9/1.5 182477*1e-9 182477*1e-9*1.5 182477*1e-9*2];
% Ir_r_z_vals = [182477*1e-9];
% h_s_vals = [-.026*2 -.026 0 .026 .026*2];
% h_s_vals = [-.026 0 .026];
% h_s_vals = [0];
h_s_vals = [-.026];
% Is_s_xy_vals = [600000*1e-9/2 600000*1e-9/1.5 600000*1e-9 600000*1e-9*1.5 600000*1e-9*2];
Is_s_xy_vals = [600000*1e-9/2 600000*1e-9 600000*1e-9*2];
% Is_s_xy_vals = [600000*1e-9];
Is_s_z_vals = [65063*1e-9/2 65063*1e-9/1.5 65063*1e-9 65063*1e-9*1.5 65063*1e-9*2];
% Is_s_z_vals = [65063*1e-9];
% h_p_vals = [-.017*2 -.017 0 .017 .017*2];
% h_p_vals = [-.017 0 .017];
% h_p_vals = [0];
h_p_vals = [.017];
% h_d_vals = [-.01*2 -.01 0 .01 .01*2]; %not default spacing because default is so freaking small
% h_d_vals = [-.01 0 .01]; %not default spacing because default is so freaking small
% h_d_vals = [0];
% h_d_vals = [-.001];
h_d_vals = [.01];
% chord_d_vals = [.5 1 2];
chord_d_vals = [1];
setup_flyer_ideal_prop_on_bottom;

h_r_averages = zeros(size(h_r_vals,2),2);
h_r_trials = zeros(size(h_r_vals,2),1);

Ir_r_x_averages = zeros(size(Ir_r_x_vals,2),2);
Ir_r_x_trials = zeros(size(Ir_r_x_vals,2),1);

Ir_r_y_averages = zeros(size(Ir_r_y_vals,2),2);
Ir_r_y_trials = zeros(size(Ir_r_y_vals,2),1);

Ir_r_z_averages = zeros(size(Ir_r_z_vals,2),2);
Ir_r_z_trials = zeros(size(Ir_r_z_vals,2),1);

h_s_averages = zeros(size(h_s_vals,2),2);
h_s_trials = zeros(size(h_s_vals,2),1);

Is_s_xy_averages = zeros(size(Is_s_xy_vals,2),2);
Is_s_xy_trials = zeros(size(Is_s_xy_vals,2),1);

Is_s_z_averages = zeros(size(Is_s_z_vals,2),2);
Is_s_z_trials = zeros(size(Is_s_z_vals,2),1);

h_p_averages = zeros(size(h_p_vals,2),2);
h_p_trials = zeros(size(h_p_vals,2),1);

h_d_averages = zeros(size(h_d_vals,2),2);
h_d_trials = zeros(size(h_d_vals,2),1);

chord_d_averages = zeros(size(chord_d_vals,2),2);
chord_d_trials = zeros(size(chord_d_vals,2),1);


trial = 1;
h_r_num = 1;
for h_r = h_r_vals
    Ir_r_x_num = 1;
    for Ir_r_x = Ir_r_x_vals
        Ir_r_y_num = 1;
        for Ir_r_y = Ir_r_y_vals
            Ir_r_z_num = 1;
            for Ir_r_z = Ir_r_z_vals
                h_s_num = 1;
                for h_s = h_s_vals
                    Is_s_xy_num = 1;
                    for Is_s_xy = Is_s_xy_vals
                        Is_s_z_num = 1;
                        for Is_s_z = Is_s_z_vals
                            h_p_num = 1;
                            for h_p = h_p_vals
                                h_d_num = 1;
                                for h_d = h_d_vals
                                    chord_d_num = 1;
                                    for chord_d_pre = chord_d_vals
%                                         chord_d = [.105 .105 .099 .094 .086 .077 .065 .05 .03]*chord_d_pre;
                                        chord_d = [.152 .152 .152 .152 .152 .152 .152 .152 .152]*chord_d_pre;
                                        
                                        Ir_r(1,1) = Ir_r_x;
                                        Ir_r(2,2) = Ir_r_y;
                                        Ir_r(3,3) = Ir_r_z;
                                        Is_s(1,1) = Is_s_xy;
                                        Is_s(2,2) = Is_s_xy;
                                        Is_s(3,3) = Is_s_z;

                                        %Calculated Parameters
                                        S_s_f = [0 0 h_s];
                                        S_r_f = [0 0 h_r];
                                        area_p = R_p*R_p*pi;
                                        area_d = R_d*R_d*pi;
                                        bladeProperties_p = 0;
    %                                     bladeProperties_p(1) = -h_p-h_r;
                                        bladeProperties_p(1) = -h_p;
                                        bladeProperties_p(2) = R_p;
                                        bladeProperties_p(3) = size(beta_p,2);
                                        bladeProperties_p = [bladeProperties_p beta_p];
                                        bladeProperties_p = [bladeProperties_p chord_p];
                                        bladeProperties_d1 = 0;
    %                                     bladeProperties_d1(1) = -h_d-h_s;
                                        bladeProperties_d1(1) = -h_d;
                                        bladeProperties_d1(2) = R_d1;
                                        bladeProperties_d1(3) = size(beta_d,2);
                                        bladeProperties_d1 = [bladeProperties_d1 beta_d];
                                        bladeProperties_d1 = [bladeProperties_d1 chord_d];
                                        bladeProperties_d2 = 0;
    %                                     bladeProperties_d2(1) = -h_d-h_s;
                                        bladeProperties_d2(1) = -h_d;
                                        bladeProperties_d2(2) = R_d2;
                                        bladeProperties_d2(3) = size(beta_d,2);
                                        bladeProperties_d2 = [bladeProperties_d2 beta_d];
                                        bladeProperties_d2 = [bladeProperties_d2 chord_d];
                                        m = m_s+m_r;
                                        I_tot = (m_s*(S_s_f')*S_s_f+m_r*(S_r_f')*S_r_f);

                                        %Run the simulation
                                        X = [3.02; 0; 0; 0; 0; 0; 0; -340; 49; .1; 0; 0; 0; 0; 0; 0; 0; -2.524];
                                        [tout Xout] = yimFlyerLite(X,6);

                                        %Analyze
                                        trials(trial,1) = max(Xout(:,10))-min(Xout(:,10));
                                        trials(trial,2) = max(Xout(:,11))-min(Xout(:,11));

                                        %populate
                                        trials(trial,3) = h_r;
                                        h_r_averages(h_r_num,1) = h_r_averages(h_r_num,1) + trials(trial,1);
                                        h_r_averages(h_r_num,2) = h_r_averages(h_r_num,2) + trials(trial,2);
                                        h_r_trials(h_r_num) = h_r_trials(h_r_num) + 1;

                                        trials(trial,4) = Ir_r_x;
                                        Ir_r_x_averages(Ir_r_x_num,1) = Ir_r_x_averages(Ir_r_x_num,1) + trials(trial,1);
                                        Ir_r_x_averages(Ir_r_x_num,2) = Ir_r_x_averages(Ir_r_x_num,2) + trials(trial,2);
                                        Ir_r_x_trials(Ir_r_x_num) = Ir_r_x_trials(Ir_r_x_num) + 1;

                                        trials(trial,5) = Ir_r_y;
                                        Ir_r_y_averages(Ir_r_y_num,1) = Ir_r_y_averages(Ir_r_y_num,1) + trials(trial,1);
                                        Ir_r_y_averages(Ir_r_y_num,2) = Ir_r_y_averages(Ir_r_y_num,2) + trials(trial,2);
                                        Ir_r_y_trials(Ir_r_y_num) = Ir_r_y_trials(Ir_r_y_num) + 1;

                                        trials(trial,6) = Ir_r_z;
                                        Ir_r_z_averages(Ir_r_z_num,1) = Ir_r_z_averages(Ir_r_z_num,1) + trials(trial,1);
                                        Ir_r_z_averages(Ir_r_z_num,2) = Ir_r_z_averages(Ir_r_z_num,2) + trials(trial,2);
                                        Ir_r_z_trials(Ir_r_z_num) = Ir_r_z_trials(Ir_r_z_num) + 1;

                                        trials(trial,7) = h_s;
                                        h_s_averages(h_s_num,1) = h_s_averages(h_s_num,1) + trials(trial,1);
                                        h_s_averages(h_s_num,2) = h_s_averages(h_s_num,2) + trials(trial,2);
                                        h_s_trials(h_s_num) = h_s_trials(h_s_num) + 1;

                                        trials(trial,8) = Is_s_xy;
                                        Is_s_xy_averages(Is_s_xy_num,1) = Is_s_xy_averages(Is_s_xy_num,1) + trials(trial,1);
                                        Is_s_xy_averages(Is_s_xy_num,2) = Is_s_xy_averages(Is_s_xy_num,2) + trials(trial,2);
                                        Is_s_xy_trials(Is_s_xy_num) = Is_s_xy_trials(Is_s_xy_num) + 1;

                                        trials(trial,9) = Is_s_z;
                                        Is_s_z_averages(Is_s_z_num,1) = Is_s_z_averages(Is_s_z_num,1) + trials(trial,1);
                                        Is_s_z_averages(Is_s_z_num,2) = Is_s_z_averages(Is_s_z_num,2) + trials(trial,2);
                                        Is_s_z_trials(Is_s_z_num) = Is_s_z_trials(Is_s_z_num) + 1;

                                        trials(trial,10) = h_p;
                                        h_p_averages(h_p_num,1) = h_p_averages(h_p_num,1) + trials(trial,1);
                                        h_p_averages(h_p_num,2) = h_p_averages(h_p_num,2) + trials(trial,2);
                                        h_p_trials(h_p_num) = h_p_trials(h_p_num) + 1;

                                        trials(trial,11) = h_d;
                                        h_d_averages(h_d_num,1) = h_d_averages(h_d_num,1) + trials(trial,1);
                                        h_d_averages(h_d_num,2) = h_d_averages(h_d_num,2) + trials(trial,2);
                                        h_d_trials(h_d_num) = h_d_trials(h_d_num) + 1;
                                        
                                        trials(trial,12) = chord_d_pre;
                                        chord_d_averages(chord_d_num,1) = chord_d_averages(chord_d_num,1) + trials(trial,1);
                                        chord_d_averages(chord_d_num,2) = chord_d_averages(chord_d_num,2) + trials(trial,2);
                                        chord_d_trials(chord_d_num) = chord_d_trials(chord_d_num) + 1;

                                        trial = trial+1
                                        chord_d_num = chord_d_num + 1;
                                    end
                                    h_d_num = h_d_num + 1
                                end
                                h_p_num = h_p_num + 1
                            end
                            Is_s_z_num = Is_s_z_num + 1
                        end
                        Is_s_xy_num = Is_s_xy_num + 1
                    end
                    h_s_num = h_s_num + 1
                end
                Ir_r_z_num = Ir_r_z_num + 1
            end
            Ir_r_y_num = Ir_r_y_num + 1
        end
        Ir_r_x_num = Ir_r_x_num + 1
    end
    h_r_num = h_r_num + 1
end

%process
h_r_averages(:,1) = h_r_averages(:,1)./h_r_trials(:);
h_r_averages(:,2) = h_r_averages(:,2)./h_r_trials(:);

Ir_r_x_averages(:,1) = Ir_r_x_averages(:,1)./Ir_r_x_trials(:);
Ir_r_x_averages(:,2) = Ir_r_x_averages(:,2)./Ir_r_x_trials(:);

Ir_r_y_averages(:,1) = Ir_r_y_averages(:,1)./Ir_r_y_trials(:);
Ir_r_y_averages(:,2) = Ir_r_y_averages(:,2)./Ir_r_y_trials(:);

Ir_r_z_averages(:,1) = Ir_r_z_averages(:,1)./Ir_r_z_trials(:);
Ir_r_z_averages(:,2) = Ir_r_z_averages(:,2)./Ir_r_z_trials(:);

h_s_averages(:,1) = h_s_averages(:,1)./h_s_trials(:);
h_s_averages(:,2) = h_s_averages(:,2)./h_s_trials(:);

Is_s_xy_averages(:,1) = Is_s_xy_averages(:,1)./Is_s_xy_trials(:);
Is_s_xy_averages(:,2) = Is_s_xy_averages(:,2)./Is_s_xy_trials(:);

Is_s_z_averages(:,1) = Is_s_z_averages(:,1)./Is_s_z_trials(:);
Is_s_z_averages(:,2) = Is_s_z_averages(:,2)./Is_s_z_trials(:);

h_p_averages(:,1) = h_p_averages(:,1)./h_p_trials(:);
h_p_averages(:,2) = h_p_averages(:,2)./h_p_trials(:);

h_d_averages(:,1) = h_d_averages(:,1)./h_d_trials(:);
h_d_averages(:,2) = h_d_averages(:,2)./h_d_trials(:);

chord_d_averages(:,1) = chord_d_averages(:,1)./chord_d_trials(:);
chord_d_averages(:,2) = chord_d_averages(:,2)./chord_d_trials(:);

%plot
figure(1);
hold off;
plot(h_r_vals,h_r_averages(:,1),'r');
hold on;
plot(h_r_vals,h_r_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('h_r');

figure(2);
hold off;
plot(Ir_r_x_vals,Ir_r_x_averages(:,1),'r');
hold on;
plot(Ir_r_x_vals,Ir_r_x_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('Ir_r_x');

figure(3);
hold off;
plot(Ir_r_y_vals,Ir_r_y_averages(:,1),'r');
hold on;
plot(Ir_r_y_vals,Ir_r_y_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('Ir_r_y');

figure(4);
hold off;
plot(Ir_r_z_vals,Ir_r_z_averages(:,1),'r');
hold on;
plot(Ir_r_z_vals,Ir_r_z_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('Ir_r_z');

figure(5);
hold off;
plot(h_s_vals,h_s_averages(:,1),'r');
hold on;
plot(h_s_vals,h_s_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('h_s');

figure(6);
hold off;
plot(Is_s_xy_vals,Is_s_xy_averages(:,1),'r');
hold on;
plot(Is_s_xy_vals,Is_s_xy_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('Is_s_xy');

figure(7);
hold off;
plot(Is_s_z_vals,Is_s_z_averages(:,1),'r');
hold on;
plot(Is_s_z_vals,Is_s_z_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('Is_s_z');

figure(8);
hold off;
plot(h_p_vals,h_p_averages(:,1),'r');
hold on;
plot(h_p_vals,h_p_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('h_p');

figure(9);
hold off;
plot(h_d_vals,h_d_averages(:,1),'r');
hold on;
plot(h_d_vals,h_d_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('h_d');

figure(10);
hold off;
plot(chord_d_vals,chord_d_averages(:,1),'r');
hold on;
plot(chord_d_vals,chord_d_averages(:,2),'g');
legend('phi pk-pk','theta pk-pk');
title('chord_d');
