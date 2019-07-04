function [Fb, Mb, aoa1] = blade(Vcg,omg,nu,pitch1,Cl, Cd, h, beta, chord, span, dSpan, rho)
    drSteps = length(beta);

    Fb = zeros(drSteps,3);
    Mb = zeros(drSteps,3);

    %calculate relative wind
    Ut1 = Vcg(1) + omg(2)*h - span*omg(3); %correct
    Up1 = Vcg(3) + span*omg(1) - nu;
    rw = atan2(Up1,Ut1);

    %lift and drag
    aoa1=mod(pi+beta+pitch1+rw,2*pi)-pi;
    ind = int32((aoa1+pi)*1000000+1);
    qs = rho/2*(Ut1.*Ut1 + Up1.*Up1).*chord;
    l1 = qs.*Cl(ind)';
    d1 = qs.*abs(Cd(ind))'; 

    %p and t components of force (n is up, t is back)
    n1 = l1.*cos(rw)+d1.*sin(rw);
    t1 = -l1.*sin(rw)+d1.*cos(rw);

    %TODO::should be done for between elements, not each element
    %force vector
    Fb(:,1) = -t1*dSpan;
    Fb(:,3) = -n1*dSpan;

    %moment vector
    Mb(:,1) = -span.*n1*dSpan;
    Mb(:,2) = -h.*t1*dSpan; %correct (sign change for h)
    Mb(:,3) = span.*t1*dSpan;
end