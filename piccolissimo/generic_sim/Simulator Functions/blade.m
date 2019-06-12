function [Fb, Mb] = blade(Vcg,omg,nu,pitch1,bladeProperties,rho)
    
    % Sorry for the globals, but loading these take too long
    global Cl_segmented Cd_segmented
    
    h = bladeProperties(1);
    R = bladeProperties(2);
    drSteps = bladeProperties(3);
    beta = bladeProperties(4:3+drSteps);
    chord = bladeProperties(4+drSteps:3+drSteps+drSteps);
    dR = R/drSteps;

    Fb = [0 0 0];
    Mb = [0 0 0];

    blade_step = 1:drSteps;
    r=dR*(blade_step-1);

    %calculate relative wind
    Ut1 = Vcg(1) + omg(2)*h - r*omg(3); %correct
    Up1 = Vcg(3) + r*omg(1) - nu;
    rw = atan2(Up1,Ut1);

    %lift and drag
    aoa1=mod(pi+beta(blade_step)+pitch1+rw,2*pi)-pi;
    ind = int32((aoa1+pi)*1000000+1);
    qs = rho/2*(Ut1.*Ut1 + Up1.*Up1).*chord(blade_step);
    l1 = qs.*Cl_segmented(ind)';
    d1 = qs.*abs(Cd_segmented(ind))'; 

    %p and t components of force (n is up, t is back)
    n1 = l1.*cos(rw)+d1.*sin(rw);
    t1 = -l1.*sin(rw)+d1.*cos(rw);

    %force vector
    Fb(1) = -sum(t1)*dR;
    Fb(3) = -sum(n1)*dR;

    %moment vector
    Mb(1) = -sum(r.*n1)*dR;
    Mb(2) = -sum(h.*t1)*dR; %correct (sign change for h)
    Mb(3) = sum(r.*t1)*dR;
end