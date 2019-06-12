function [F_p, M_p, F_d, M_d] = ComputeAero()
    % Computes aerodynamic forces in the flyer frame

    % Load flyer properties
    global B_d B_p Rb_f Rr_s Rr_f Rs_b angle_s angle_r Vcg omg omg_b omg_r nu pitch1_d pitch_p bladeProperties_d1 bladeProperties_p rho
    % drag blades computed in body frame
    F_d = [0 0 0];
    M_d = [0 0 0];
    for drag_blade = 1:B_d
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*drag_blade/B_d; 
        Rb_f = [cos(angle_s + blade_angle), -sin(angle_s + blade_angle), 0; sin(angle_s + blade_angle), cos(angle_s + blade_angle), 0; 0, 0, 1];
        [F, M] = blade(Rb_f*Vcg,Rb_f*(omg+[0; 0; omg_b]),nu,pitch1_d,bladeProperties_d1,rho);
        F_d = F_d + (F*Rb_f);
        M_d = M_d + (M*Rb_f);
    end

    % propeller blades computed in rotor frame
    F_p = [0 0 0];
    M_p = [0 0 0];
    for prop_blade = 1:B_p
        % end with 2*pi so future calculations don't care about the blade angle
        blade_angle = 2*pi*prop_blade/B_p;
        Rr_s = [cos(angle_r + blade_angle), -sin(angle_r + blade_angle), 0; sin(angle_r + blade_angle), cos(angle_r + blade_angle), 0; 0, 0, 1];
        Rr_f = Rr_s*Rs_b*Rb_f;
%         [F, M] = blade(Rr_f*(Vcg),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),bladeProperties_p,rho);
        [F, M] = blade(Rr_f*(Vcg.*[0;0;1]),Rr_f*(omg+[0; 0; omg_r]),nu,pitch_p(prop_blade),bladeProperties_p,rho); %Vcg.*[0;0;1] for shrowded prop
        F_p = F_p + (F*Rr_f);
        M_p = M_p + (M*Rr_f);
    end
end