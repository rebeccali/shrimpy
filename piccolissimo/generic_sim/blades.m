function [F1 M1 F2 M2] = blades(Vcg,omg,nu,pitch1,pitch2,bladeProperties,rho)
h = bladeProperties(1);
R = bladeProperties(2);
drSteps = bladeProperties(3);
beta = bladeProperties(4:3+drSteps);
chord = bladeProperties(4+drSteps:3+drSteps+drSteps);
persistent Cd Cl
if isempty(Cd)
    load('Cd');
    load('Cl');
end

F1 = [0 0 0];
M1 = [0 0 0];
F2 = [0 0 0];
M2 = [0 0 0];

for i = 1:drSteps
    r=R*i/drSteps;
    %Blade 1
        %calculate relative wind
        Ut1 = Vcg(1) - (omg(2)*h+r*omg(3));
        Up1 = Vcg(3) + r*omg(1) - nu;
        
        %lift and drag
        aoa1=mod(pi+beta(i)+pitch1+atan2(Up1,Ut1),2*pi)-pi;
        l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa1));
        d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa1));

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
    
    %Blade 2
        %calculate relative wind
        Ut2 = -Vcg(1) + omg(2)*h-r*omg(3);
        Up2 = Vcg(3) - r*omg(1) - nu;
        
        %lift and drag
        aoa2=mod(pi+beta(i)+pitch2+atan2(Up2,Ut2),2*pi)-pi;
        l2 = rho/2*(Ut2*Ut2 + Up2*Up2)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa2));
        d2 = rho/2*(Ut2*Ut2 + Up2*Up2)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa2));

        %p and t components of force (n is up, t is back)
        n2 = l2*cos(atan2(Up2,Ut2))+d2*sin(atan2(Up2,Ut2));
        t2 = -l2*sin(atan2(Up2,Ut2))+d2*cos(atan2(Up2,Ut2));

        %force vector
        F2(1) = F2(1)+t2*(R/drSteps);
        F2(3) = F2(3)-n2*(R/drSteps);
        
        %moment vector
        M2(1) = M2(1) + (r)*n2*(R/drSteps);
        M2(2) = M2(2) - h*t2*(R/drSteps);
        M2(3) = M2(3) + (r)*t2*(R/drSteps);
end

