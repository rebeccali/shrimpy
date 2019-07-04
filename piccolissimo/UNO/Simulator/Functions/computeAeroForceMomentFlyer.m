function [forcesAeroBody_f, momentsAeroBody_f, forcesProp_f, momentsProp_f aoa_body, aoa_prop] = computeAeroForceMomentFlyer()
%% Compute aerodynamic forces in the flyer frame
%  Returns the forces in the flyer frame, moments in the flyer frame,
%  angle of attack of the body blades, angle of attack of prop blades

% B_p is the number of propellor blades (shaft propellor)
% B_d is the number of body drag plates

    %% First Calculate Body drag plates
    % Load flyer properties
    % drag blades computed in body frame
    beta_d =
    angle_s = X(14);


    F_d = zeros(length(beta_d),3);
    M_d = zeros(length(beta_d),3);
    aoa_body = zeros(length(beta_d),bodyParams.B_d);

    for drag_blade = 1:bodyParams.B_d
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*drag_blade/bodyParams.B_d;
        Rb_f = [cos(angle_s + blade_angle), -sin(angle_s + blade_angle), 0; sin(angle_s + blade_angle), cos(angle_s + blade_angle), 0; 0, 0, 1];

        [F, M, aoa] = blade(Rb_f*Vcg,Rb_f*(omg+[0; 0; omg_b]),nu,pitch1_d,Cl_drag,Cd_drag,h_d,beta_d,chord_d,span_d,dSpan_d,rho);

        F_d = F_d + F*Rb_f;
        M_d = M_d + M*Rb_f;
        aoa_body(:,drag_blade) = aoa;
    end

    %% Calculate Propellor drag plates
    % propeller blades computed in rotor frame
    F_p = zeros(length(beta_p),3);
    M_p = zeros(length(beta_p),3);
    aoa_prop = zeros(length(beta_p),propParams.B_p);
    for prop_blade = 1:propParams.B_p
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*prop_blade/propParams.B_p;
        Rr_s = [cos(angle_r + blade_angle), -sin(angle_r + blade_angle), 0; sin(angle_r + blade_angle), cos(angle_r + blade_angle), 0; 0, 0, 1];
        Rr_f = Rr_s*Rs_b*Rb_f;
%         [F, M, aoa] = blade(Rr_f*(Vcg.*[0;0;1]),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),Cl_prop,Cd_prop,h_p,beta_p,chord_p,span_p,dSpan_p,rho); %Vcg.*[0;0;1] for shrowded prop
        [F, M, aoa] = blade(Rr_f*(Vcg),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),Cl_prop,Cd_prop,h_p,beta_p,chord_p,span_p,dSpan_p,rho); % no shrowded prop

        F_p = F_p + F*Rr_f;
        M_p = M_p + M*Rr_f;
        aoa_prop(:,prop_blade) = aoa;
    end

end
