function [F_p, M_p, F_d, M_d, aoa_p, aoa_d] = ComputeAero()
    % Computes aerodynamic forces in the flyer frame

    % Load flyer properties
    global B_d B_p Rb_f Rr_s Rr_f Rs_b angle_s angle_r Vcg omg omg_b omg_r nu pitch1_d pitch_p h_d beta_d chord_d span_d dSpan_d h_p beta_p chord_p span_p dSpan_p rho Cl_prop Cd_prop Cl_drag Cd_drag
    % drag blades computed in body frame
    F_d = zeros(length(beta_d),3);
    M_d = zeros(length(beta_d),3);
    aoa_d = zeros(length(beta_d),B_d);
    for drag_blade = 1:B_d
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*drag_blade/B_d; 
        Rb_f = [cos(angle_s + blade_angle), -sin(angle_s + blade_angle), 0; sin(angle_s + blade_angle), cos(angle_s + blade_angle), 0; 0, 0, 1];
        [F, M, aoa] = blade(Rb_f*Vcg,Rb_f*(omg+[0; 0; omg_b]),nu,pitch1_d,Cl_drag,Cd_drag,h_d,beta_d,chord_d,span_d,dSpan_d,rho);
        
        F_d = F_d + F*Rb_f;
        M_d = M_d + M*Rb_f;
        aoa_d(:,drag_blade) = aoa;
    end

    % propeller blades computed in rotor frame
    F_p = zeros(length(beta_p),3);
    M_p = zeros(length(beta_p),3);
    aoa_p = zeros(length(beta_p),B_p);
    for prop_blade = 1:B_p
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*prop_blade/B_p;
        Rr_s = [cos(angle_r + blade_angle), -sin(angle_r + blade_angle), 0; sin(angle_r + blade_angle), cos(angle_r + blade_angle), 0; 0, 0, 1];
        Rr_f = Rr_s*Rs_b*Rb_f;
%         [F, M, aoa] = blade(Rr_f*(Vcg.*[0;0;1]),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),Cl_prop,Cd_prop,h_p,beta_p,chord_p,span_p,dSpan_p,rho); %Vcg.*[0;0;1] for shrowded prop
        [F, M, aoa] = blade(Rr_f*(Vcg),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),Cl_prop,Cd_prop,h_p,beta_p,chord_p,span_p,dSpan_p,rho); % no shrowded prop

        F_p = F_p + F*Rr_f;
        M_p = M_p + M*Rr_f;
        aoa_p(:,prop_blade) = aoa;
    end
end