function [F1, M1] = blade(Vcg,omg,nu,pitch1,bladeProperties,rho)
h = bladeProperties(1);
R = bladeProperties(2);
drSteps = bladeProperties(3);
beta = bladeProperties(4:3+drSteps);
chord = bladeProperties(4+drSteps:3+drSteps+drSteps);

persistent Cd_segmented Cl_segmented
if isempty(Cd_segmented)
    load('Cd_segmented');
    load('Cl_segmented');
end

F1 = [0 0 0];
M1 = [0 0 0];

for i = 1:drSteps
    r=R*(i-1)/drSteps;
    %calculate relative wind
%     Ut1 = Vcg(1) - (omg(2)*h+r*omg(3)); %old
    Ut1 = Vcg(1) + omg(2)*h - r*omg(3); %correct
    Up1 = Vcg(3) + r*omg(1) - nu;

    %lift and drag
    aoa1=mod(pi+beta(i)+pitch1+atan2(Up1,Ut1),2*pi)-pi;
    l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(Cl_segmented(int32((aoa1+pi)*100+1)));
    d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(Cd_segmented(int32((aoa1+pi)*100+1)));
%     l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa1));
%     d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa1));

    %p and t components of force (n is up, t is back)
    n1 = l1*cos(atan2(Up1,Ut1))+d1*sin(atan2(Up1,Ut1));
    t1 = -l1*sin(atan2(Up1,Ut1))+d1*cos(atan2(Up1,Ut1));

    %force vector
    F1(1) = F1(1)-t1*(R/drSteps);
    F1(3) = F1(3)-n1*(R/drSteps);

    %moment vector
    M1(1) = M1(1) - r*n1*(R/drSteps);
    M1(2) = M1(2) + h*t1*(R/drSteps);
    M1(3) = M1(3) + r*t1*(R/drSteps);
end

