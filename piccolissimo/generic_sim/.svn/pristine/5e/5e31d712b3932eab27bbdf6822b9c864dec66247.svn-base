%Test script
%Runs many simulations and evals the difference
    trial = 0;
    global  drSteps h_r m_r Ir_r h_s m_s Is_s h_p R_p beta_p chord_p H_p h_d R_d beta_d chord_d H_d rho g S_s_f S_r_f area_p area_d bladeProperties_p bladeProperties_d1 bladeProperties_d2 m I_tot counter

%Changable Parameters
        disp('Setting parameters');
    %Changable Parameters
        %Rotor properties
        h_r = -.068; %Height of rotor cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_r = .063; %rotor mass in kg
        Ir_r = [203536.211032,87.172143,0.003051;87.172143,26959.742608,-0.127636;0.003051,-0.127636,182477.134918]*1e-9; %Rotor inertia at rotor cg in flyer frame
%         Ir_r = [203536.211032,0,0;0,26959.742608,0;0,0,182477.134918]*1e-9; %Rotor inertia at rotor cg in flyer frame

        %Stator properties
        h_s = .026; %Height of stator cg below flyer cg (cg assumed axial) (meters) (positive is down)
        m_s = .164; %stator mass in kg
        Is_s = [624513.130008,70.172180,379.489186;70.172180,587292.296260,625.057297;379.489186,625.057297,65063.603653]*1e-9; %Stator inertia at stator cg in flyer frame
%         Is_s = [624513.130008,0,0;0,587292.296260,0;0,0,65063.603653]*1e-9; %Stator inertia at stator cg in flyer frame

        %Propeller properties
        h_p = -.017; %Height of propeller below ROTOR cg (meters) (positive is down)
%         h_p = -h_r; %Height of propeller below ROTOR cg (meters) (positive is down)
        R_p = .178; %Single blade radius (meters)
%         beta_p = 0.2*ones(1,10); %Propeller twist (relative to zero lift)(assumed constant) (radians)
%         beta_p = [0 0 0.0499 0.1810 0.1788 0.1474 0.1184 0.0913 .062]; %printed
        beta_p = [0 0 .397 .341 .2998 .258 .224 .194 .171 .151 .136]; %14x6
%         chord_p = .02*ones(1,10); %Chord lengh (assumed constant) (meters)
%         chord_p =[.01 .01 .015 .023 .023 .021 .019 .017 .015]; %printed
        chord_p =[.0228 .0171 .00247 .0269 .0299 .0317 .0313 .0291 .0249 .0205 .0146]; %14x6
        H_p = R_p; %approx height above prop that air is moved (estimated to be radius of prop)

        %Drag plate properties
%         h_d = .001; %Height of center of drag plate below stator cg (meters) (positive is down)
        h_d = -h_s; %Height of center of drag plate below stator cg (meters) (positive is down)
        R_d1 = .140; %Single plate radius (meters)
        R_d2 = .140; %Single plate radius (meters)
%         beta_d = pi/2*ones(1,drSteps); %Propeller twist (relative to zero lift)(assumed constant) (radians)
        beta_d = [pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %to be function of r
%         chord_d = .15*ones(1,drSteps); %Chord lengh (assumed constant) (meters)
        chord_d =[.105 .105 .099 .094 .086 .077 .065 .05 .03]; %to be function of r
        H_d = 0; %approx height above prop that air is moved (estimated to be zero)

        %Environment
        rho = 1.225; %air density (sea level 1.225) (kg/m^3)
        g = 9.8; %gravity acceleration (m/s^2)
        
        close(figure(9));
        
%         h_d_vals = -.078:.013:.078;
%     for h_d = h_d_vals
    h_p_vals = -.136*2:.017:.136*2;
    for h_p = h_p_vals
        trial = trial+1;

    %Calculated Parameters
        S_s_f = [0 0 h_s];
        S_r_f = [0 0 h_r];
        area_p = R_p*R_p*pi;
        area_d = R_d*R_d*pi;
        bladeProperties_p = 0;
        bladeProperties_p(1) = -h_p-h_r;
        bladeProperties_p(2) = R_p;
        bladeProperties_p(3) = size(beta_p,2);
        bladeProperties_p = [bladeProperties_p beta_p];
        bladeProperties_p = [bladeProperties_p chord_p];
        bladeProperties_d1 = 0;
        bladeProperties_d1(1) = -h_d-h_s;
        bladeProperties_d1(2) = R_d1;
        bladeProperties_d1(3) = size(beta_d,2);
        bladeProperties_d1 = [bladeProperties_d1 beta_d];
        bladeProperties_d1 = [bladeProperties_d1 chord_d];
        bladeProperties_d2 = 0;
        bladeProperties_d2(1) = -h_d-h_s;
        bladeProperties_d2(2) = R_d2;
        bladeProperties_d2(3) = size(beta_d,2);
        bladeProperties_d2 = [bladeProperties_d2 beta_d];
        bladeProperties_d2 = [bladeProperties_d2 chord_d];
        m = m_s+m_r;
        I_tot = (Is_s+Ir_r+m_s*(S_s_f')*S_s_f+m_r*(S_r_f')*S_r_f);
        counter = 0;
        
    %Run the simulation
        X = [3; 0; 0; 0; 0; 0; 0; -340; 72; .1; 0; 0; 0; 0; 0; 0; 0];
        [tout Xout] = yimFlyerLite(X,15);
        
    %Analyze
        dt = diff(tout);
        dt(size(diff(tout),1)+1) = 0;
        av(trial,1) = sum(abs(Xout(:,10)'*dt));
        av(trial,2) = sum(abs(Xout(:,11)'*dt));
        peak(trial,1) = max(Xout(:,10))-min(Xout(:,10));
        peak(trial,2) = max(Xout(:,11))-min(Xout(:,11));
        if trial == 1
            Xout10 = Xout(:,10);
            Xout11 = Xout(:,11);
            tout_all = tout;
        else
            Xout10 = [Xout10; Xout(:,10)];
            Xout11 = [Xout11; Xout(:,11)];
            tout_all = [tout_all; tout];
        end
        figure(9);
        hold on
        plot3(tout,Xout(:,10),(h_p+h_r)*ones(size(tout)),'r');    
        plot3(tout,Xout(:,11),(h_p+h_r)*ones(size(tout)),'g'); 
        pause(1);
    end
    figure(7)
    hold off
    plot(h_p_vals+h_r,av(1:trial,:));
    hold on
    plot(h_p_vals+h_r,av(1:trial,1)+av(1:trial,2),'c');
    legend('phi','theta','sum');
    
    figure(8)
    hold off
    plot(h_p_vals+h_r,peak(1:trial,:));
    hold on
    legend('phi','theta');

    figure(9);
    view(3);
%     hold off
%     plot3(tout_all(1),Xout10(:,1),h_d_vals(1)+h_s,'r');
%     plot3(tout_all(1),Xout11(:,1),h_d_vals(1)+h_s,'g');
%     hold on
%     for i = 2:trial
%         i
%         plot3(tout_all(i),Xout10(:,i),h_d_vals(i)+h_s,'r');
%         plot3(tout_all(i),Xout11(:,i),h_d_vals(i)+h_s,'g');
%     end