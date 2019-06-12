function [F1 M1] = blade2(Vcg,omg,nu,pitch1,bladeProperties,rho)
h = bladeProperties(1);
R = bladeProperties(2);
drSteps = bladeProperties(3);
beta = bladeProperties(4:3+drSteps);
chord = bladeProperties(4+drSteps:3+drSteps+drSteps);
% persistent Cd Cl
% if isempty(Cd)
%     load('Cd');
%     load('Cl');
% end
persistent Cd_segmented Cl_segmented
if isempty(Cd_segmented)
    load('Cd_segmented');
    load('Cl_segmented');
end

F1 = [0 0 0];
M1 = [0 0 0];

R_step = R/drSteps;

for i = 1:drSteps
    r=(i-1)*R_step;
    
    %calculate relative wind
    Ut1 = Vcg(1) - (omg(2)*h+r*omg(3));
    Up1 = Vcg(3) + r*omg(1) - nu;

    %lift and drag
    rw_aoa = atan2(Up1,Ut1);
    aoa1=mod(pi+beta(i)+pitch1+rw_aoa,2*pi)-pi;
    l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(Cl_segmented(int32((aoa1+pi)*100+1)));
    d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(Cd_segmented(int32((aoa1+pi)*100+1)));
%     l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa1));
%     d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa1));

    %p and t components of force (n is up, t is back)
    c_rw_aoa = cos(rw_aoa);
    s_rw_aoa = sin(rw_aoa);
    n1 = l1*c_rw_aoa+d1*s_rw_aoa;
    t1 = -l1*s_rw_aoa+d1*c_rw_aoa;

    %force vector
    F1(1) = F1(1)-t1*(R_step);
    F1(3) = F1(3)-n1*(R_step);

    %moment vector
    M1(1) = M1(1) - r*n1*(R_step);
    M1(2) = M1(2) + h*t1*(R_step);
    M1(3) = M1(3) + r*t1*(R_step);
end

